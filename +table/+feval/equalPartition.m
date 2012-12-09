function I = equalPartition(nGroups,varargin) 
% I = table.par.equalPartition(5,A) assigns a set index to each row of A in such a
% way that:
%          i.   there are nGroups,
%          ii.  each set has roughly the same number of elements 

[K,U] = table.indexing.createKey(varargin{:});
nKeys = length(U);
f = accumarray(K,1,[nKeys,1]);
cf = cumsum(f);
avg = length(K)./nGroups;
assgn = ceil(cf./avg);
I = assgn(K);

end




