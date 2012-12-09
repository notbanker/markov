function [laggedDx,I] = differences(X,lags,varargin)
% Create a matrix of differences after grouping.
mL = max(lags)+1;
[laggedX,I] = table.laggedValues(X,(1:mL),varargin{:});
allLaggedDx = -diff([X,laggedX],1,2);
laggedDx = allLaggedDx(:,lags);
end
