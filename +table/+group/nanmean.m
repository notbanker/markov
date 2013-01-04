function avgX = nanmean(X,varargin)
f = @(x) nanmean(x);
avgX = table.group.fun(X,f,varargin{:});
end