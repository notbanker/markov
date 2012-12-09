function tf = isfield(ds,field)
% function tf = isfield(ds,field)

if isa(ds,'struct'),
    tf = isfield(ds,field);
elseif isa(ds,'dataset'),
    tf = ismember(field,ds.Properties.VarNames);
end