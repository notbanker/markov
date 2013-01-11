function [laggedX,I] = values(X,lags,varargin)
% Create a matrix of lagged values after grouping
if size(X,2)>1,
    error('Anticipating column vector');
end
[keys,uKeys] = table.indexing.createKey(varargin{:});
J = size(uKeys,1);
I = table.indexing.previousOccurances(keys,lags,J);
laggedX = nan(size(X,1),length(lags));
laggedX(~isnan(I)) = X(I(~isnan(I)));
end