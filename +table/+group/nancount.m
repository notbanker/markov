function avgX = nancount(X,varargin)
f = @(x) sum(~isnan((x)));
avgX = table.group.fun(X,f,varargin{:});
end