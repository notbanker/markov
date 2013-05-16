function [w,v,x] = hlmtrain(y,Q,varargin)

lambda = 100; % Ratio of measurement var to dynamic std

lattice = hlmlattice(y,varargin{:});
v = y(5:end)-y(1:end-4);
w = diff(y);
nIter = 20;
vAll = zeros(0,1);
wAll = zeros(0,1);

for iter=1:nIter,
    [x,p] = hlmdecode(y,w,v,lattice);
    
    % errors
    vKer = fitdist([y-x';x'-y],'kernel');
    v = [random(vKer,5000,1);trnd(1.7*ones(200,1))];
    
    wKer = fitdist([diff(x');-diff(x')],'normal');
    w = [random(wKer,5000,1);trnd(1.7*ones(10,1))];
    
    if 0,
        % balance variance
        v12 = var(w)+var(v);
        c1 = 1/(1 + lambda);
        c2 = 1/(1 + 1/lambda);
        std1 = sqrt(c1*v12);
        std2 = sqrt(c2*v12);
        w = w*std1/std(w);
        v = v*std1/std(v);
        
        vs(iter) = v12;
        
    else
       
        w = w-mean(w);
        v = v-mean(v);
        
        wExcess = std(w)/sqrt(Q);
        
        wOld = w;
        vOld = v;
        w = w/wExcess;
        
        % Increase the variance of v to make up for w's reduction
        wVarDecrease = var(wOld)-var(w);
        
        if abs(wVarDecrease)<var(v),
           multiplier = sqrt(wVarDecrease/var(v)+1);
           v = v*multiplier;
           
           shouldBeZero = var(v)+var(w)-var(vOld)-var(wOld),
        end
        ws(iter) = wExcess;
    end
    
    
    % plot
    clf; subplot(3,2,1); plot(x);hold on; plot(y,'g*');
    subplot(3,2,3); hist(w,400);xlabel('w (dynamics)');xs= axis;axis([-5,5,xs(3),xs(4)]);
    subplot(3,2,4); hist(v,400);xlabel('v (measurement)');xs= axis;axis([-5,5,xs(3),xs(4)]);
    subplot(3,2,2); plot([1:iter],ws(1:iter));
    subplot(3,2,5); normplot(w);
    subplot(3,2,6); normplot(v);
   
    drawnow;pause(0.3);
    
end
