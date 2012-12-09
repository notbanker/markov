function [K,U,R,O] = createKey(varargin)
% Create unique contiguous keys into a collection of variables
% (treated as columns in a database table) all of whom have length m
%
% K        m x 1 - Unique keys
% U(i,:)   The i'th unique row-combinations of columnwise id's keys
% R{i}     Mapping from columwise id's into unique representatives of X(:,i)
% O(k)     The position of the last occurance of key k

n = length(varargin);
m = size(varargin{1},1);
U = nan(m,n);
R = cell(n,1);
for k=1:n,
    if size(varargin{k},1)==m,
        [R{k},~,U(:,k)] = unique(varargin{k});
    else
        error('Expecting variables of the same length');
    end
end
[U,O,K] = unique(U,'rows');

end