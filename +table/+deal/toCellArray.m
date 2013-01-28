function C = toCellArray(ds,fieldsToInclude)
% Deal out tabular struct or dataset fields into cell array
assert(isa(ds,'dataset') || isa(ds,'struct'),'Expecting table as first argument');
assert(iscell(fieldsToInclude));
fns = table.fieldnames(ds);
[~,obsLoc] = ismember(fieldsToInclude,fns);
if all(obsLoc>0),
    nInclude = length(fieldsToInclude);
    C = cell(nInclude,1);
    for k=1:nInclude,
        x_ = ds.(fns{obsLoc(k)});
        C{k} = x_;
    end
else
    warning(['table.deal.toCellArray was passed an invalid field']);
    f = find(obsLoc<1);
    fieldsToInclude{f(:)},
    error('Given up');
end

end