function uds = unique(ds,varargin)
% Supply a list of fields. Returns a table where every row contains a unique combination
% of values taken by these fields (the last occurance is chosen). 
[K,U,R,O] = table.indexing.createKey(table.deal.toVarargout(ds,varargin{:}));
uds = table.select.rows(ds,O); 
end

