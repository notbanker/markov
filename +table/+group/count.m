function avgX = count(X,varargin)
f = @(x) numel(x);
avgX = table.group.fun(X,f,varargin{:});
end