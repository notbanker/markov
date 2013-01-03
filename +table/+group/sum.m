function avgX = sum(X,varargin)
f = @(x) sum(x);
avgX = table.group.fun(X,f,varargin{:});
end