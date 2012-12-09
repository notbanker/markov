function varargout = toVarargout(ds,varargin)
% Deal out dataset or struct fields into stand alone variables
fns = fieldnames(ds);
[~,obsLoc] = ismember(varargin(:),fns);
if all(obsLoc>0),
    n = length(varargin);
    varargout = cell(n,1);
    for k=1:n,
        varargout{k} = ds.(fns{obsLoc(k)});
    end
else
    warning(['Invalid field name(s) requested']);
    f = find(obsLoc<1);
    varargin{f(:)},
    error('giving up');
end
end