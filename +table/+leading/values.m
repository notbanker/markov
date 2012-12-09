function [leadingX,I] = values(X,lags,varargin)
% Create a matrix of leading values after grouping
if size(X,2)>1,
    error('Anticipating column vector');
end
[keys,uKeys] = table.indexing.createKey(varargin{:});
n = size(keys,1);
I = table.indexing.previousOccurances(keys(n:-1:1),lags,uKeys);
I = I(n:-1:1,:);
Xr = X(n:-1:1);
leadingX = nan(size(Xr,1),length(lags));
leadingX(~isnan(I)) = Xr(I(~isnan(I)));
end