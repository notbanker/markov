function I = previousOccurances(key,lags,keyMax,C)
% I=previousOccurances(key,lags) returns a matrix with length(key) rows and
% length(lags) columns. 
%   If lags(k)=1, column I(i,k) reports the previous position j satisfying key(j)==key(i). 
%   If lags(k)=2 I(i,k) reports the time before that, and so forth. 
%
% NaN's are returned if there was no previous occurance. Here is an example with lags=2,
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
% lags can be any vector of desired lags such as [1,2,4,8]
%
% The remaining two optional arguments can assist performance. 
%    keyMax must be greater than equal to the largest key value, if provided   
%    C defaults to (1..keyMax). 
%
% C can be used to report the previous position for a *different* key rather than the key which appears.
% This is mildly confusing and documented in createKeyChoice.tex. The notation 'C' is shared by 
% TABLE.LAGGED.VALUESOFANOTHER and TABLE.INDEXING.CREATEKEYCHOICE, incidentally. For an
% example of this more advanced usage see TABLE.TEST.UNITTEST_LAGGEDVALUESOFANOTHER

% TODO: This function is a bottleneck and might be rewritten as mex


if ~isvector(lags) || any(lags<=0) || any(abs(ceil(lags)-lags)>1e-6),
    error('Expecting positive integer lags');
end
[a,b] = size(key);
if b~=1,
    error('Expecting vector key');
end

if nargin>=3 && length(keyMax)>1,
    warning('You may be using an old form of this function where U was supplied instead of keyMax');
    keyMax = size(keyMax,1);
end
if nargin<3,
    keyMax = max(key); % Reason apparent in loop below
end

if nargin<4,
    C = (1:keyMax)';
end

maxLags = max(lags);
prevI = nan(keyMax,maxLags);
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
    warning('Using a slower loop. It is strongly recommended that the third parameter (keyMax) be supplied.');
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
