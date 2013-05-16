classdef roll
    
    properties (Constant)
       worm = true;
       nSteps = 100;
       nLattice = 600;
       nu = 0.8;  % Bid-offer distance
       v = 2.1;   % Degrees of freedom for t-distribution
       R = 2.3^2; % Measurement error variance
       Q = 0.25^2; % Dynamics 
       kalMultiplier = 8; % Amount to increase Q by in the Kalman filter equations
       percentOfMassInFlatTail = 0.1;
       minimumEmitProbability = 1e-8; 
       nLags = 5;      % How far to look back
       nLook = 3;  % How much to cheat (look forward)
       plotDance = false; % Show the evolving posterior state distribution
       floorGainAtZero = false;
       kappa = 0.002;   % Mean reversion in process (helps stay on grid)
       addSurprisesToSim = true;
       addJumpToMeasurements = true;
       showFakeData = false;
       showHead = false;
    end
    
    % Analysis of the roll model and estimators thereof
    
    methods (Static) % Roll model simulation
        
        function [q,x,y] = rollSim(n,Q,R,v,nu,kappa,x0)
          % Q  variance per time step
          % R  variance in measurement error
          % nu jump in measurement
            
            if nargin<3,
                R = 1;
            end
            if nargin<2,
                dxStd = 1;
            else
                dxStd = sqrt(Q);
            end
            
            x = nan(n,1);
            x(1) = dxStd*trnd(v)+x0;
            for k=2:n,
                x(k) = x(k-1) -kappa*x(k-1) + dxStd*randn(1);
                if n>100 && k==100,
                    x(50) = x(50)+5;
                end
                if n>80 && k==80,
                    x(80) = x(80)-5;
                end
                if roll.addSurprisesToSim && rand(1)<0.05,
                  x(k) = x(k) + trnd(1.2); 
                end
               
            end
            % More jumps for fun
            
            q = ceil(rand(n,1)*2)*2-3;
            [~,tvar] = tstat(v);
            noise = trnd(v*ones(n,1))*sqrt(R)/tvar;
            if roll.addJumpToMeasurements,
               for k=1:length(noise),
                  if rand(1)<0.3,
                     noise(k) = noise(k)+10*rand(1)-5; 
                  end
               end
            end
            y = x+q*nu+noise;
            
        end
        
    end
    
    methods (Static)  % Estimation & filtering
        
        
        function results = decode(y,Q,R,v,nu,m)
            
            % Infer state from univariate observations
            %
            % y    - observation vector
            % Q    - variance (dynamics) 
            % R    - variance (measurements)
            % nu   - bid/offer gap
            % m    - roughly half the lattice size
            %
            % Returns a struct with fields:
            %
            % results.p  - posterior probabilities
            % results.xPlus - posterior mean
            % results.plusP - posterior variance
            
            maxR = 3*R;
            minR = R/3;
            
            %% Define filtering transition model
            n = length(y);  
            x0 = 1.5;
           
            %% Choose a model 
            %[both,grid] = discretizer.coinAndGauss(x0,Q,R,nu,m,n);
            [both,grid] = discretizer.coinAndT(x0,Q,R,v,nu,m,n);
            enlargement = 1;
            while min(grid)>min(y)-2*sqrt(R) || max(grid)<max(y)+2*sqrt(R),
               enlargement = enlargement*1.2,
               [both,grid] = discretizer.coinAndT(x0,Q*enlargement,R*enlargement,v,nu*enlargement,m,n*enlargement);
            end
        
            % Add extreme jumps
            both = both+roll.percentOfMassInFlatTail/m;
            both = both/sum(both);
           
    
            %% Checks
            results.discreteE = sum(grid.*both);
            results.discreteR = sum((grid-results.discreteE).^2.*both);
            s0 = m+1; % the middle 
            em = roll.emit(both,s0);
            drift = discretizer.atomicGauss(grid,x0,sqrt(Q));
            Edrift = sum(grid.*drift);
            results.discreteQ = sum((grid-Edrift).^2.*drift);
            s0 = m+1;
            tr = roll.emit(drift,s0);
            [results.yNearest,seq] = discretizer.nearest(y,grid);
            
            %% Filter
            [results.p,results.pWorm] = roll.hmmposterior(seq',tr,em);
            results.xPlus = grid*results.p;
            results.plusP = grid.^2*results.p - results.xPlus.^2;
            results.yBefore = [y(1);y(2:end)];
            results.xMinus = [results.xPlus(1),results.xPlus(1:end-1)];
            results.xBefore = [results.xMinus(1),results.xMinus(1:end-1)];
            for sampleNo=1:n,
                pw = squeeze(results.pWorm(:,:,sampleNo)); 
                for t=1:roll.nLook,
                   results.xWorm(:,sampleNo) = grid*pw;
                end
            end
          
            
            if roll.plotDance,
            figure;
            for k=1:size(results.p,2),plot(results.p(:,k)); drawnow;pause(0.3);end
            end
            
            %% Create Fake observations to fool regular Kalman filter
            results.yFake = y;
            for k=3:n,
                xMinus = results.xPlus(k-1);
                xPlus = results.xPlus(k);
                minusP = results.plusP(k-1)+Q;
                plusP = results.plusP(k);
                impliedK = 1-plusP/minusP;
                if impliedK<0.01,
                    impliedR = maxR;
                elseif impliedK>1,
                    impliedR = minR;
                else
                    impliedR = (1/impliedK-1)*minusP;
                end
                if impliedR>maxR,
                    impliedR = maxR;
                elseif impliedR<minR,
                    impliedR = minR;
                end
                if roll.floorGainAtZero && impliedK<0,
                    K = 0;
                    results.yImplied(k) = y(k);
                else
                    K = minusP/(minusP+impliedR);
                    results.yImplied(k) = (xPlus-(1-K)*xMinus)/K;
                end
                results.impliedR(k) = impliedR;
                results.impliedK(k) = impliedK;
                results.K(k) = K;
            end
            
            
            
             
        end
       
        function [p,pWorm] = hmmposterior(seq,tr,em)
            n = size(tr,1);
            nSamples = length(seq);
            p = zeros(n,nSamples);
            pWorm = zeros(n,roll.nLook,nSamples);
            for sampleNo=1:nSamples-roll.nLook,
               start = max(sampleNo-roll.nLags,1);
               lah = min(roll.nLook,nSamples-sampleNo);
               seq_ = seq(start:sampleNo+lah);
               p_ = hmmdecode(seq_,tr,em);
               p(:,sampleNo) = p_(:,end-lah);
               pWorm(:,:,sampleNo) = p_(:,end-lah+1:end);
            end
            
        end
        
        
        function [xHat,yHat,xPlus,yPlus] = kalman(y,Q,R,x0)
            % Kalman filter for 
            x = x0;
            P = 0.0001;
            n = length(y);
            H = 1;
            if isscalar(R),
               R = R*ones(n,1); 
            end
            for k=1:n
                xHat(k) = x;
                yHat(k) = H*x;
                P_ = P+Q;
                x_ = x;
                K = P_*H'/(H*P_*H'+R(k));
                if K>1 || K<0,
                   error('WTF?!'); 
                end
                x  = x_ + K*(y(k)-H*x_);
                P  = (1-K*H)*P_;
                xPlus(k) = x;
                yPlus(k) = H*x;
               
            end
            
            if 0,
                plot(1:n,y,'r*');hold on;
                plot(1:n,xPlus);hold on;
            end
        end
        
    end
    
    
    
    
    methods (Static) % Emission and translation matrices
        
        
        function em = emit(em0,s0)
            % Generate translation invariant hidden Markov model
            n = length(em0);
            m = (n-1)/2;
            if ~m==ceil(m),
                error('Expecting odd number of emissions');
            end
            em = zeros(n,n);
            for s1=1:n,
                em(s1,:) = shiftEmit(em0,s0,s1)';
            end
            
            
            em = em+roll.minimumEmitProbability;
            em = discretizer.normalizeRows(em);
            
            function em1 = shiftEmit(em0,s0,s1)
                % em0(j) is the emission probability to signal j given that we are in state s0
                % em1(j) is the emission probability to signal j given that we are in state s1
                if s1==s0,
                    em1 = em0;
                elseif s0>s1,
                    n = length(em0);
                    em1 = zeros(n,1);
                    em1(1:n-(s0-s1)) = em0(1+s0-s1:n);
                elseif s1>s0,
                    n = length(em0);
                    em1 = zeros(n,1);
                    em1(s1-s0+1:end) = em0(1:n-(s1-s0));
                end
            end
        end
        
        
    end
    
    methods(Static) % Unit tests
      
     
       
       function unitTest_hlmtrain
           % Illustrates the use of observation adjustments that allow the
           % Kalman filter to act like a non-linear filter
           
           close all;
           
           %% Simulate roll model
           n = roll.nSteps;
           m = roll.nLattice;
           Q = roll.Q;
           R = roll.R;
           x0 = 1.3;
           nu = roll.nu;
           v = roll.v; 
           kappa = roll.kappa; 
           
           
           if v<2.1,
              error('Things won''t work too well if measurement error distribution has no second moment ?!'); 
           end
           
           [~,x,y] = roll.rollSim(n,Q,R,v,nu,kappa,x0);
           
           dx = diff(x);
           estQ = var(dx(dx<4*std(dx)));
           [w,v,x] = hlmtrain(y,Q);
       end
       
         
       function unitTest_obsManager
           % Illustrates the use of observation adjustments that allow the
           % Kalman filter to act like a non-linear filter
           
           close all;
           
           %% Simulate roll model
           n = roll.nSteps;
           m = roll.nLattice;
           Q = roll.Q;
           R = roll.R;
           x0 = 1.3;
           nu = roll.nu;
           v = roll.v; 
           kappa = roll.kappa; 
           
           
           if v<2.1,
              error('Things won''t work too well if measurement error distribution has no second moment ?!'); 
           end
           
           [~,x,y] = roll.rollSim(n,Q,R,v,nu,kappa,x0);
           
           %% Filter numerically, creating fake observations
           res = roll.decode(y,Q,R,v,nu,m);
           xHmmPlus = res.xPlus;
           x0 = xHmmPlus(1);
           
           %% Filter real observations (Kalman)
           [xMinus,~,xPlus] = roll.kalman(y,roll.kalMultiplier*res.discreteQ,res.discreteR,x0);
            
           %% Filter fake observations (Kalman)
           [zMinus,~,zPlus] = roll.kalman(res.yImplied,res.discreteQ,res.impliedR,x0);
           
           
           
           %% Plots
           close all;
           figure;
           for sampleNo=3:n-roll.nLook,
               zW = res.xWorm(:,sampleNo);
               clf;
               replot(xPlus,xHmmPlus,zPlus,zW,sampleNo);drawnow;pause(0.1);
               xlabel('Time');
               ylabel('Price');
               legend('Data','Kalman Filter','Human');
           end
           
           
           %% Throw out burn-in and report errors
           xHmmMinus = [xMinus(1),xHmmPlus(1:end-1)];
           nBurn = ceil(0.2*n);
           x_ = x(nBurn:end-1);
           xMinus_ = xMinus(nBurn:end-1);
           xPlus_ = xPlus(nBurn:end-1);
           zMinus_ = zMinus(nBurn:end-1);
           xHmmMinus_ = xHmmMinus(nBurn:end-1);
           xHmmPlus_ = xHmmPlus(nBurn:end-1);
           y_ = y(nBurn:end-1);
           zPlus_ = zPlus(nBurn:end-1);
           impliedK_ = res.impliedK(nBurn:end-1);
           
           rms = @(x) num2str(ceil(100*sqrt(mean(x.^2)))/100);
           dist = @(x) num2str(ceil(100*mean(abs(diff(x))))/100);
           bias = @(x) num2str(ceil(100*mean(x))/100);
           
           xlabel(['K/H/F RMS=',rms(x_-xMinus_'),'/',rms(x_-xHmmMinus_'),'/',rms(x_-zMinus_'),'  ',...
                  'D=',dist(xMinus_),'/',dist(xHmmMinus_),'/',dist(zMinus_),' ',...
                  ' BIAS=',bias(xMinus_'-x_),'/',bias(xHmmMinus_'-x_),'/',bias(zMinus_'-x_)]);
              
           %% Take a look at the fake observations   
           innov = y_-xHmmMinus_';
           move = xHmmPlus_'-xHmmMinus_';
           
           
           figure;
           subplot(2,2,1);
           plot(innov,move,'.');xlabel('Innovation');ylabel('State adjustment');grid;
           subplot(2,2,2);
           plot(innov,impliedK_','.');xlabel('Innovation');ylabel('K');grid;
           
          
           function  replot(xK,xH,zP,zW,sampleNo)
               xK = xK(1:sampleNo);
               xH = xH(1:sampleNo);
               zP = zP(1:sampleNo);
               
               n = length(xK);
               ahead = (n+(1:roll.nLook))';
               plot(1:length(y),y,'b*');hold on;grid;
               plot(1:n,xK,'c-','LineWidth',3,'MarkerSize',3);hold on;
               plot(1:n,xH,'m','LineWidth',3,'MarkerSize',3);hold on;
               if roll.showHead,
                   plot(ahead',zW','m-','LineWidth',2,'MarkerSize',2);hold on;
               end
               if roll.showFakeData,
                   plot(1:n,res.yImplied,'ro');
               end
           end
           
              
       end
          
       
      
          
   end
   
   
  
    
end