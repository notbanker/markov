function [laggedDy,I] = differencesOfAnother(X,lags,f,varargin)
% Create a matrix of lagged differences after grouping, but reporting the lagged difference
% of choice(j) rather than the lagged value of j
mL = max(lags)+1;
[laggedY,I] = table.lagged.valuesOfAnother(X,(1:mL),f,varargin{:});
allLaggedDy = -diff([X,laggedY],1,2);
laggedDy = allLaggedDy(:,lags);
end
