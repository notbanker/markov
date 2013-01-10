function uds = unique(ds,fieldsToGroupBy,fieldsToShow)
% Supply a list of fields. Returns a table where every row contains a unique combination
% of values taken by these fields (the last occurance is chosen). 
vals = table.deal.toNumericArray(ds,fieldsToGroupBy);

[K,U,R,O] = table.indexing.createKey(  );
uds = table.select.rows(ds,O); 
if nargin>=3,
   uds = table.select.cols(uds,fieldsToShow); 
end
end

