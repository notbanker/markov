function I = equalPartition(nSets,varargin) 
% I = equalPartition(5,A) assigns a set index to each row of A in such a
% way that:
%          i.   there are nSets,
%          ii.  each set has roughly the same number of elements 

[K,U] = obsLagged.createKey(varargin{:});
nKeys = length(U);
f = accumarray(K,1,[nKeys,1]);
cf = cumsum(f);
avg = length(K)./nSets;
assgn = ceil(cf./avg);
I = assgn(K);

end




