function fX = groupFun(X,f,varargin)
% Apply function to all subgroups as defined by varargin
% The function f should take vectors to scalars or vectors to a vector of the same length
keys_ = table.createKey(varargin{:});
[uKeys,~,keys] = unique(keys_);
fX = nan(size(X));
for k=1:length(uKeys),
    I = (keys==uKeys(k));
    x = feval(f,X(I));
    fX(I) = x;
end
end
