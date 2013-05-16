classdef discretizer
    
    % Utilities useful in simple univariate filtering problems where dynamics and measures
    % are represented on a grid of evenly spaced lattice points. See
    % hmldecode for example
    
    methods (Static) % Miscellaneous crud
       
       function y = normalize(x)
           if isvector(x),
               y = x/nansum(x);
           else
               error('Cannot normalize matrix, etc');
           end
       end
       
       function y = normalizeRows(x)
          y = nan(size(x));
          for k=1:size(x,1),
             y(k,:) = discretizer.normalize(x(k,:)); 
          end
       end
        
   end
    
    methods(Static) % Discretization, nearest neighbour, categorization, quantiles
        
        function [yNear,I] = nearest(y,lattice)
            % Return position of nearest lattice point
            n = length(y);
            m = length(lattice);
            v = 1:m;
            I = interp1(lattice,v,y,'nearest','extrap');
            yNear = lattice(I);
        end
        
        function [zDis,f,g] = discretize(z,unit)
            % Discretizes and returns function handles for the maps
            nBelow = ceil((-min(z)+unit/2)/unit);
            nAbove = ceil((max(z)-unit/2)/unit);
            edges = (-nBelow:nAbove)*unit+unit/2;
            if min(z)<min(edges),
                error('Urgh');
            end
            if max(z)>max(edges),
                error('Urgh');
            end
            edgesT = edges(:);
            f = @(x) discretizer.categorize(z,edges);
            g = @(n) edgesT(n)-unit/2;
            zDis = f(z);
            if any(zDis==0),
                error('Groan');
            end
        end
        
        function zDis = categorize(z,edges)
            [m,b] = size(z);
            if b==1,
                E = [-inf,edges(:)'];
                EDGES = ones(m,1)*E;
                nE = length(E);
                Z = z*ones(1,nE);
                zDis = sum(Z>=EDGES,2);
            else
                error('Z should be a vector');
            end
            
        end
        
        
        function f = quantileRepresentative(x,qntl)
            x1 = quantile(x,min(1,1.05*qntl));
            x2 = quantile(x,0.95*qntl);
            f = find(x<=x1 & x>=x2,1,'first');
            if isempty(f),
                error('Could not find percentile');
            end
        end
        
        
    end
    
    methods (Static) % Elementary distributions
        
        function w = atomicGauss(lattice,mu,sig)
            f = @(x) normpdf(x,mu,sig);
            w = discretizer.atomic(lattice,f);
        end
        
        function w = atomicT(lattice,mu,sig,v)
            f = @(x) tpdf((x-mu)/sig,v);
            w = discretizer.atomic(lattice,f);
        end
        
        
        function w = coin(lattice,mu,q)
            w = 0.5*discretizer.atom(lattice,mu+q)+0.5*discretizer.atom(lattice,mu-q);
        end
       
        
        function w = atomic(lattice,f)
            % Find weights best approximating any distribution when mass must be distributed at known lattice points
            if isa(f,'function_handle'),
                nFine = 10000;
                x = linspace(lattice(1),lattice(end),nFine);
                fx = f(x);
                dx = x(2)-x(1);
                fx = fx/sum(fx);
                w = zeros(size(lattice));
                for k=1:nFine
                    wx = discretizer.atom(lattice,x(k));
                    w = w+fx(k)*wx;
                end
            elseif isa(f,'double');
                % Samples
                x = f;
                w = zeros(size(lattice));
                nSamples = length(x);
                for k=1:nSamples,
                    wx = discretizer.atom(lattice,x(k));
                    w = w+wx/nSamples;
                end
                
            end
        end
        
        
        function w = atom(lattice,x)
            % Find weights best approximating an atomic distribution centred at
            % x, when mass must be distributed at known lattice points
            n = length(lattice);
            w = zeros(size(lattice));
            [yep,loc] = ismember(x,lattice);
            if yep,
                w(loc(1))=1;
            elseif all(lattice<x),
                w(1) = 1;
            elseif all(lattice>x),
                w(n) =1;
            else
                fHigh = find(lattice>=x,1,'first');
                fLow = find(lattice<=x,1,'last');
                xHigh = lattice(fHigh(1));
                xLow = lattice(fLow(1));
                aHigh = (x-xLow)./(xHigh-xLow);
                aLow = 1-aHigh;
                w(fLow) = aLow;
                w(fHigh) = aHigh;
            end
            
        end
        
    end
    
    methods (Static) % Projection, Convolution
        
        function w = project(wOld,oldGrid,lattice)
            % Move from one representation of a measure to another
            w = zeros(size(lattice));
            for k=1:length(wOld),
               xk = oldGrid(k);
               wk_ = wOld(k)*discretizer.atom(lattice,xk);
               w = w + wk_;
            end
            
            oldE = sum(wOld.*oldGrid);
            newE = sum(w.*lattice);
            bias = newE-oldE;
            if abs(bias)>0.0001,
               warning(['Projection introduced bias of ',num2str(bias)]);
            end
            
            
        end
         
        function [w,lattice,wAB,latticeAB] = convolve(wA,latticeA,wB,latticeB,doProjectOnB)
           
            % Given two atomic measures representing variables A and B,
            % produce an atomic measure wAB representing A+B on new, larger
            % latticeAB and then project back to lattice=latticeA (by default) or latticeB
            
            
            %% Ensure odd and compatible
            nA = length(wB);
            nB = length(wA);
            mA = (nA-1)/2;
            mB = (nB-1)/2;
            if ceil(mA)~=mA || ceil(mB)~=mB,
               error('Anticipating odd length atomic measures'); 
            end
            dxB = latticeB(2)-latticeB(1);
            if any(abs(diff(latticeB)-dxB)>0.00001), 
                error('Anticipating equaly spaced latticeB');
            end
            dxA = latticeA(2)-latticeA(1);
            if any(abs(diff(latticeA)-dxA)>0.00001), 
                error('Anticipating equaly spaced latticeA');
            end
            if abs(dxA-dxB)>0.0001,
                error('Anticipating compatible lattices');
            end
            a = latticeA(1)-dxA;
            b = latticeB(1)-dxB; 
            
            wA_ = [wA(1),wA(1:end-1)];  % So that wA_(i-j+1) = wA(i-j)
            wAB = conv(wB,wA_);         % wAB(k) = sum  wB(j) wA_(k-j+1) 
                                        %        = sum  wB(j) wA(k-j)
                                        %        = P[B=b+j*dx] P(A = a+ (k-j)*dx)
                                        %        = P[A+B = a + b + k*dx]
            latticeAB = (1:length(wAB))*dxA+a+b;
            
            % Check bias
            EA = sum(wA.*latticeA);
            EB = sum(wB.*latticeB);
            EAB = sum(wAB.*latticeAB);
            bias = EAB-EA-EB;
            if abs(bias)>0.001,
               warning(['Convolution introduced bias of ',num2str(bias)]); 
            end
            
            % Project
            if nargin>=5 && doProjectOnB,
                lattice = latticeB;
            else
                lattice = latticeA;
            end
            w = discretizer.project(wAB,latticeAB,lattice);
            
        end
        
    end
    
    methods (Static) % Particular combinations of distributions
        
         function [both,lattice] = coinAndGauss(x0,Q,R,nu,m,n)
            lattice = discretizer.setLatticeFromVariance(x0,nu.^2+3*R+n*Q,m);
            if nu>0.01,
                jump = discretizer.coin(lattice,0,nu);  % Mean zero jump
                meas = discretizer.atomicGauss(lattice,x0,sqrt(R));
                both = discretizer.convolve(jump,lattice,meas,lattice);
            else
                both = discretizer.atomicGauss(lattice,x0,sqrt(R));
            end
            
            if 0,
                plot(lattice,jump,lattice,meas,lattice,both); legend('jump','meas','both');
            end
        end
        
        function [both,lattice] = coinAndT(x0,Q,R,v,nu,m,n)
            
            lattice = discretizer.setLatticeFromVariance(x0,2*nu.^2+3*R+1.2*n*Q,m);
            if nu>0.01,
                jump = discretizer.coin(lattice,0,nu);  % Mean zero jump
                meas = discretizer.atomicT(lattice,x0,sqrt(R),v);
                both = discretizer.convolve(jump,lattice,meas,lattice);
                if 1,
                     plot(lattice,jump,lattice,meas,lattice,both); legend('jump','meas','both');
                end
            else
                both = discretizer.atomicT(lattice,x0,sqrt(R),v);
            end
            
            
        end
        
         function lattice = setLatticeFromVariance(x1,R,m)
           % Create lattice points based on initial value and variance
           devos = 3*(-m:m)/m;
           lattice = x1+sqrt(R)*devos;
        end
        
        
    end
      
    methods (Static) % Unit tests
                
       function okay = unitTest_atom
           grid = linspace(1,134,17);
           x = 34/3;
           em = roll.atom(grid,x);
           okay = abs(x-sum(em.*grid))<0.0001;
       end
       
       
       function okay = unitTest_discretize
           try
               x = randn(100,1);
               unit = 0.5;
               [xDis,f,g] = discretizer.discretize(x,unit);
               xDisCheck = f(xDis);
               xCheck = g(xDis(:));
               okay1 = all(abs(xCheck-x)<=unit/2);
               okay2 = all(xDisCheck==xDis);
               okay = okay1 && okay2;
           catch
               okay = false;
           end
       end
        
        function okay = unitTest_atomicGauss
            lattice = linspace(-6,10,100);
            w = discretizer.atomicGauss(lattice,1.3,2);
            okay = abs(sum(lattice.*w)-1.3)<0.01;
        end
        
        
        
    end
    
end