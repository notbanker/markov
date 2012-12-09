function [laggedY,I] = valuesOfAnother(X,lags,f,varargin)
% Create a matrix of lagged values after grouping, but reporting the lagged value
% of choice(j) rather than the lagged value of j
[K,U,R,C] = table.indexing.createKeyChoice(f,varargin{:});
I = table.indexing.previousOccurances(K,lags,U,C);
laggedY = nan(size(X,1),length(lags));
laggedY(~isnan(I)) = X(I(~isnan(I)));
end