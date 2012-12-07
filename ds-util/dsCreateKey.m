function key = dsCreateKey(ds,varargin)
% Create unique contiguous key given field names
vars = cell(length(varargin),1);
[vars{:}] = dsDeal(ds,varargin{:});
key = createKey(vars{:});
end


