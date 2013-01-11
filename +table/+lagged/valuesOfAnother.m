function [laggedY,I] = valuesOfAnother(X,lags,f,varargin)
% function [laggedY,I] = valuesOfAnother(X,lags,f,varargin)
% Create a matrix of lagged values after grouping, but reporting the lagged value
% of choice(j) rather than the lagged value of j
if isa(f,'function_handle') || isa(f,'char'),
    [K,U,~,C] = table.indexing.createKeyChoice(f,varargin{:});
    I = table.indexing.previousOccurances(K,lags,size(U,1),C);
    laggedY = nan(size(X,1),length(lags));
    laggedY(~isnan(I)) = X(I(~isnan(I)));
else
    error('You may have forgotten the f argument');
end
end