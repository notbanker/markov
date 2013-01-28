function [leadingY,I] = valuesOfAnother(X,lags,f,varargin)
% Create a matrix of leading values after grouping, but reporting the leading value
% of choice(j) rather than the leading value of j
[K,U,R,C] = table.indexing.createKeyChoice(f,varargin{:});
n = length(K);
nUniqueKeys = size(U,1);
I = table.indexing.previousOccurances(K(n:-1:1),lags,nUniqueKeys,C);
I = I(n:-1:1,:);
Xr = X(n:-1:1);
leadingY = nan(size(Xr,1),length(lags));
leadingY(~isnan(I)) = Xr(I(~isnan(I)));
end