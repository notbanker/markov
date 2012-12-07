classdef obsLagged
    
    % Creation of lagged variables by groups
    
    properties(Constant)
       trainingWheelsOn = false;  
    end
    
    methods (Static) % 'Table' manipulation functions operating on variables
        
        function avgX = groupAverage(X,varargin)
            f = @(x) nanmean(x);
            avgX = obsLagged.groupFun(X,f,varargin{:});
        end
        
        function avgX = groupCount(X,varargin)
            f = @(x) sum(~isnan(x));
            avgX = obsLagged.groupFun(X,f,varargin{:});
        end
        
        function avgX = groupSum(X,varargin)
            f = @(x) nansum(x);
            avgX = obsLagged.groupFun(X,f,varargin{:});
        end
        
        function fX = groupFun(X,f,varargin)
            % Apply function to all subgroups as defined by varargin
            % The function f should take vectors to scalars or vectors to a vector of the same length
            keys_ = obsLagged.createKey(varargin{:});
            [uKeys,~,keys] = unique(keys_);
            fX = nan(size(X));
            for k=1:length(uKeys),
                I = (keys==uKeys(k));
                x = feval(f,X(I));
                fX(I) = x;
            end
        end
        
        function [laggedDx,I] = laggedDifferences(X,lags,varargin)
            % Create a matrix of differences after grouping. 
            mL = max(lags)+1;
            [laggedX,I] = obsLagged.laggedValues(X,(1:mL),varargin{:});
            allLaggedDx = -diff([X,laggedX],1,2);
            laggedDx = allLaggedDx(:,lags);
        end

        function [laggedX,I] = laggedValues(X,lags,varargin)
            % Create a matrix of lagged values after grouping 
            if size(X,2)>1,
               error('Anticipating column vector'); 
            end
            [keys,uKeys] = obsLagged.createKey(varargin{:});
            I = obsLagged.previousOccurances(keys,lags,uKeys);
            laggedX = nan(size(X,1),length(lags));
            laggedX(~isnan(I)) = X(I(~isnan(I)));
        end
        
        function [leadingX,I] = leadingValues(X,lags,varargin)
            % Create a matrix of leading values after grouping 
            if size(X,2)>1,
               error('Anticipating column vector'); 
            end
            [keys,uKeys] = obsLagged.createKey(varargin{:});
            n = size(keys,1);
            I = obsLagged.previousOccurances(keys(n:-1:1),lags,uKeys);
            I = I(n:-1:1,:);
            Xr = X(n:-1:1);
            leadingX = nan(size(Xr,1),length(lags));
            leadingX(~isnan(I)) = Xr(I(~isnan(I)));
        end

        function [laggedDy,I] = laggedDifferencesOfAnother(X,lags,f,varargin)
            % Create a matrix of lagged differences after grouping, but reporting the lagged difference
            % of choice(j) rather than the lagged value of j
            mL = max(lags)+1;
            [laggedY,I] = obsLagged.laggedValuesOfAnother(X,(1:mL),f,varargin{:});
            allLaggedDy = -diff([X,laggedY],1,2);
            laggedDy = allLaggedDy(:,lags);
        end
        
        function [laggedY,I] = laggedValuesOfAnother(X,lags,f,varargin)
            % Create a matrix of lagged values after grouping, but reporting the lagged value
            % of choice(j) rather than the lagged value of j
            [K,U,R,C] = obsLagged.createKeyChoice(f,varargin{:});
            I = obsLagged.previousOccurances(K,lags,U,C);
            laggedY = nan(size(X,1),length(lags));
            laggedY(~isnan(I)) = X(I(~isnan(I)));
        end
        
        function [leadingY,I] = leadingValuesOfAnother(X,lags,f,varargin)
            % Create a matrix of leading values after grouping, but reporting the leading value
            % of choice(j) rather than the leading value of j
            [K,U,R,C] = obsLagged.createKeyChoice(f,varargin{:});
            n = length(K);
            I = obsLagged.previousOccurances(K(n:-1:1),lags,U,C);
            I = I(n:-1:1,:);
            Xr = X(n:-1:1);
            leadingY = nan(size(Xr,1),length(lags));
            leadingY(~isnan(I)) = Xr(I(~isnan(I)));
        end
        
        function I = nextOccurances(key,leads,varargin)
           n = size(key,1);
           reversedKey = key(n:-1:1);
           reversedI = obsLagged.previousOccurances(reversedKey,leads,varargin{:});
           I = reversedI(n:-1:1);
        end

        function I = previousOccurances(key,lags,U,C)
            % I=previousOccurances(key,lags) for vectors key and lags 
            % returns a matrix with length(key) rows and 
            % length(lags) columns. If lags(k)=1, column I(i,k) reports the
            % previous position j satisfying key(j)==key(i).If lags(k)=2 I(i,k)
            % reports the time before that, and so forth. NaN's are
            % returned if there was no previous occurance. 
            %
            % Example with lags=2, 
            %
            %      key   I  
            %      2   NaN
            %      4   NaN
            %      4   NaN
            %      1   NaN
            %      2   NaN
            %      2     1
            %      3   NaN
            %      3   NaN
            %      4     2
            %      2     5
            %      3     7
            %      3     8
            %
            % lags can be any vector of desired lags such as [1,2,4,8] etc
            %
            % C is optional and defaults to (1..tildeJ) where tildeJ is the number of 
            % unique keys. It can be used to report the previous position for a *different* key rather than
            % the key which appears. For example, we might want to report
            % the last value of a customer buy each time a dealer trade
            % occurs. See createKeyChoice and accompanying latex
            % documentation. 
            
            % TODO: This loop is the bottleneck and might be rewritten as mex
            J = size(U,1);
            if nargin<4,
               C = (1:J)'; 
            end
            
            if ~isvector(lags) || any(lags<=0) || any(abs(ceil(lags)-lags)>1e-6),
                error('Expecting positive integer lags');
            end
            [a,b] = size(key);
            if b~=1,
                error('Expecting vector key');
            end
            if nargin>=3,
                nKinds = length(U);
            else
                nKinds = length(unique(key(1:min(500,a))))+10;
            end
            maxLags = max(lags);
            prevI = nan(nKinds,maxLags);
            I = nan(a,length(lags));
            if 0, % alternative algo, slightly faster because vectorized, but still clobbers cache?
                  % ** doesn't handle C() yet **
                uniqKeys = unique(sort(key));
                for k=uniqKeys'
                    w = key==k;  % <-- The bottleneck.
                    f = find(w);
                    for i=1:length(lags)
                        L = lags(i);
                        I(f(1+L:end),i) = f(1:end-L);
                    end
                end
                return;
                Inew = I;
            end
            if nargin>=3,
                for k=1:a,
                    I(k,:) = prevI(C(key(k)),lags);
                    l_ = prevI(key(k),1:(maxLags-1));
                    prevI(key(k),2:maxLags) = l_;
                    prevI(key(k),1) = k;
                end
            else
                warning('Third parameter is highly recommended for speed');
                % Loop is repeated here because it might be slower
                for k=1:a,
                    if key(k)>size(prevI,1),
                        prevI = [prevI;nan(key(k)+10-size(prevI,1),maxLags)]; % expand lookup if necessary
                    end
                    I(k,:) = prevI(C(key(k)),lags);
                    % roll down
                    l_ = prevI(key(k),1:(maxLags-1));
                    prevI(key(k),2:maxLags) = l_;
                    prevI(key(k),1) = k;
                end
            end
%             agree= isequalwithequalnans(I,Inew);
%             if agree, disp(['yay lags=' num2str(lags)]);
%             else
%                 disp(['boo lags=' num2str(lags)]); 
%                 error('stopping for debug'); 
%             end
        end
        
    end
      
    methods (Static) % Indexing/key creation utilities
     
        function [K,U,R,O] = createKey(varargin)
            % Create unique contiguous keys into a collection of variables
            % (treated as columns in a database table) all of whom have length m
            %
            % K        m x 1 - Unique keys 
            % U(i,:)   The i'th unique row-combinations of columnwise id's keys 
            % R{i}     Mapping from columwise id's into unique representatives of X(:,i)
            % O(k)     The position of the last occurance of key k 
          
            n = length(varargin);
            m = size(varargin{1},1);
            U = nan(m,n);
            R = cell(n,1);
            for k=1:n,
                if size(varargin{k},1)==m,
                    [R{k},~,U(:,k)] = unique(varargin{k});
                else
                    error('Expecting variables of the same length');
                end
            end
            [U,O,K] = unique(U,'rows');
            
        end
        
        function [tildeK,tildeU,tildeR,phiF] = createKeyChoice(f,varargin)
            % Create unique keys into a collection of equal length
            % variables (treated as table columns), and lift the user
            % supplied mapping f to the level of these keys. 
            % 
            % For a detailed explanation see obsLagged.tex
            % 
            % 
            % tildeK      - Unique key for rows of X 
            % tildeU(i,:) - The i'th unique row-combinations of column primary keys 
            % tildeR{i}   - Mapping from columwise primary keys into unique representatives of X(:,i) and f(X)(:,i)
            % phiF(k)     - The key of f(X(i,:)) when X(i,:) has key k 
            
            % Step 1: Create unique key (i.e. lookup) into rows of X
            [K,U,R,O] = obsLagged.createKey(varargin{:});
            J = size(U,1);
            
            % Step 2: Apply user choice function f and re-index
            n = length(varargin);
            allVars = cell(n,1);  
            [allVars{:}] = applyAndAppend(f,U,R,O,varargin{:});
            [hatK,tildeU,tildeR,hatO] = obsLagged.createKey(allVars{:});
            tildeK = hatK(K);
            
            if obsLagged.trainingWheelsOn,
                % Verify that old and new indexes agree
                m = length(varargin{1});
                for m_=1:m,
                    for n_=1:n,
                        if isnumeric(tildeR{n_}) || islogical(tildeR{n_}),
                            vTilde =  tildeR{n_}(tildeU(tildeK(m_),n_));
                            v =  R{n_}(U(K(m_),n_));
                            assert(isequal(v,vTilde));
                        else
                            vTilde =  tildeR{n_}{tildeU(tildeK(m_),n_)};
                            v =  R{n_}{U(K(m_),n_)};
                            assert(isequal(v,vTilde));
                        end
                    end
                end
                % Verify definitions
                obsLagged.assertRightInverse(K,O,(1:length(K)));
                obsLagged.assertRightInverse(hatK,hatO,(1:length(hatK)));
            end
            
            % Step 3: Lift the user choice function to the new keys
            F = [((J+1):2*J)';(1:J)'];             
            invHatK = obsLagged.rightInverse(hatK,1:J);
            phiF = obsLagged.compose(hatK,F,invHatK);
            
            if obsLagged.trainingWheelsOn,
                assert(isequalwithequalnans(hatK(F(1:J)),phiF(hatK(1:J))));
                invHatKCheck = obsLagged.rightInverseSimple(hatK(1:J));
                assert(isequalwithequalnans(invHatKCheck,invHatK));
                obsLagged.assertRightInverse(hatK(1:J),invHatK,(1:J));
            end
            
            
            function varargout = applyAndAppend(f,U,R,O,varargin)
                
                % Preallocate varargout
                n = length(varargin);
                nJ = length(O);
                nUniqueKeys = length(U);
                for k=1:n,
                    if isnumeric(varargin{k}) || islogical(varargin{k}),
                        varargout{k} = [varargin{k}(O);zeros(nJ,1)];
                    else
                        varargout{k} = [varargin{k}(O);cell(nJ,1)];
                    end
                end
                
                % In the second half, overwrite the value for the first variable
                for j=1:nUniqueKeys,
                    vars = cell(n,1);
                    for vNo=1:n,
                        if isnumeric(R{vNo}) || islogical(R{vNo}),
                            vars{vNo} = R{vNo}(U(j,vNo));
                        else
                            vars{vNo} = R{vNo}{U(j,vNo)};
                        end
                    end
                    varsout_ = cell(n,1);
                    [varsout_{:}] = feval(f,vars{:});
                    for k=1:n,
                        if isnumeric(varargout{k}) || islogical(varargout{k}),
                            varargout{k}(nJ+j) = varsout_{k};
                        else
                            varargout{k}{nJ+j} = varsout_{k};
                        end
                    end
                end
            end
        end
        
        function Z = compose(varargin)
            switch nargin,
                case 0,
                    Z = [];
                case 1,
                    Z = varargin{1};
                case 2,
                    Z = composeTwo(varargin{1},varargin{2});
                otherwise
                    tail = composeTwo(varargin{end-1},varargin{end});
                    Z = obsLagged.compose(varargin{1:end-2},tail);
            end
            
            function XY = composeTwo(X,Y)
                XY = nan(size(Y));
                legitY = ~isnan(Y) & Y>0;
                XY(legitY) = X(Y(legitY));
            end
        
        end
        
        function R = rightInverseSimple(L)
           R = nan(nanmax(L),1);
           for j=1:length(R),
              f = find(L==j,1,'first');
              if isempty(f),
                  R(j) = nan;
              else
                  R(j) = f(1); 
              end
           end
            
        end
        
        function R = rightInverse(L,dom)
           rangeL = L(dom);
           R = nan(nanmax(rangeL),1);
           isInDom = ismember((1:length(L))',dom);
           for j=1:length(R),
              f = find(L==j & isInDom,1,'first');
              if isempty(f),
                  R(j) = nan;
              else
                  R(j) = f(1); 
              end
           end
            
        end
        
        function assertRightInverse(L,R,dom)
            rangeL = L(dom);
            assert(~any(isnan(R(rangeL))));
            for y = rangeL(:)',
                assert(L(R(y))==y);
            end
        end
      
    end
    
    methods(Static) % Testing
        
        function ds = dsExample1
            ds = dataset;
            ds.name = {'sam','bill','mary'}';
            ds.age = [15,13,12]';
            ds.gender = {'male','male','female'}';
        end
        
        function ds = dsExample2(nSamples)
            ds = dataset;
            ds.temperature = ceil(50+rand(nSamples,1)*20);
            ds.dayOfWeek = rem((1:nSamples)',7);
            ds.rain = ceil(rand(nSamples,1)*2);
            ds.dayOfWeekName = cell(nSamples,1);
            for k=1:nSamples,
               ds.dayOfWeekName{k} = datestr(ds.dayOfWeek(k),'D');
            end
        end
           
        function [dowName1,rain1] = monday(dowName,rain)
            switch lower(dowName)
                case 's'
                    dowName1 = 'F';
                    rain1 = rain+1;
                otherwise
                    dowName1 = 'T';
                    rain1 = rain-1;
            end
        end
        
        function [report] = accumulate(X,varargin)
            % Add X grouping on other variables and ignoring nans
            [key,U,R] = obsLagged.createKey(varargin{:});
            X(isnan(X))=0;
            valid = double(~isnan(X));
            report.sum = accumarray(key,X,[length(U),1]);
            sumXX = accumarray(key,X.^2,[length(U),1]);
            report.count = accumarray(key,valid,[length(U),1]);
            report.mean = report.sum./report.count;
            report.var = sumXX./report.count;
            report.std = sqrt(report.var);
            report.table = obsLagged.reps2table(report.std,U,R);
        end
           
        function report = reps2table(X,U,R)
            % Generate a report based on accumulate output
            m = size(X,1);
            n = length(R);
            [mCheck,nCheck] = size(U);
            assert(nCheck==n && mCheck==m,'Invalid dimensions');
            report = cell(m,n+1);
            for n1=1:n,
                if iscell(R{n1}),
                    for m1=1:m,
                        report{m1,n1} = R{n1}{U(m1,n1)};
                    end
                else
                    for m1=1:m,
                        report{m1,n1} = R{n1}(U(m1,n1));
                    end
                end
            end
            for m1=1:m,
               report{m1,n1+1} = X(m1); 
            end
        end
            
        
        function unitTest_createKeyChoice
            ds = obsLagged.dsExample2(7);
            f = @obsLagged.monday; 
            [key,U,R,C] = obsLagged.createKeyChoice(f,ds.dayOfWeekName,ds.rain);
        end
        
        function unitTest_laggedValuesOfAnother
            ds = obsLagged.dsExample2(15);
            lags = (1:2)';
            ds.laggedTemperature = obsLagged.laggedValues(ds.temperature,lags,ds.dayOfWeekName);
            
            f = @obsLagged.onMonTueWedReportTue;
            [ds.strangeLaggedTemperature,I] = obsLagged.laggedValuesOfAnother(ds.temperature,lags,f,ds.dayOfWeekName);
            
            g = @obsLagged.onMonTueWedReportTue_Rain; 
            [ds.rainLaggedTemperature,I] = obsLagged.laggedValuesOfAnother(ds.temperature,lags,g,ds.dayOfWeekName,ds.rain);
         
            ds,
        end
        
        function unitTest_leadingValuesOfAnother
            ds = obsLagged.dsExample2(15);
            lags = (1:2)';
            ds.leadingTemperature = obsLagged.leadingValues(ds.temperature,lags,ds.dayOfWeekName);
            
            f = @obsLagged.onMonTueWedReportTue;
            [ds.strangeLeadingTemperature,I] = obsLagged.leadingValuesOfAnother(ds.temperature,lags,f,ds.dayOfWeekName);
            
            g = @obsLagged.onMonTueWedReportTue_Rain; 
            [ds.rainLeadingTemperature,I] = obsLagged.leadingValuesOfAnother(ds.temperature,lags,g,ds.dayOfWeekName,ds.rain);
         
            ds,
        end
        
        
        function [dowName1] = onMonTueWedReportTue(dowName)
            % An example of a user choice function
            switch dowName
                case {'M','T','W'},
                    dowName1 = 'T';
                otherwise
                    dowName1 = dowName;
            end
        end
        
        function [dowName1,rain] = onMonTueWedReportTue_Rain(dowName,rain)
            % An example of a user choice function
            rain1 = rain;
            switch dowName
                case {'M','T','W'},
                    dowName1 = 'T';
                otherwise
                    dowName1 = dowName;
            end
        end
        
        
    end
 
    
end