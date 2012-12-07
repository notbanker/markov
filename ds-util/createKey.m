function K = createKey(varargin)
% When supplied a list of variables interpreted as columns in a database,
% K=CREATEKEY(A,B,C) returns a unique key
n = length(varargin);
m = size(varargin{1},1);
for k=1:n,
    if size(varargin{k},1)==m,
        [~,~,U(:,k)] = unique(varargin{k});
    else
        error('Expecting variables of the same length');
    end
end
[~,~,K] = unique(U,'rows');
end