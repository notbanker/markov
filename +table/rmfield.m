function ds = rmfield(ds,field)
% function ds = rmfield(ds,field)

if isa(ds,'struct'),
    ds = rmfield(ds,field);
elseif isa(ds,'dataset'),
    if ischar(field) && ismember(field,ds.Properties,VarNames),
        ds.(field) = [];
    elseif iscell(field),
        isValid = ismember(field,ds.Properties.VarNames);
        field = field(isValid);
        for k=1:length(field),
            ds.(field{k}) = [];
        end
    end
else
    error('Expecting struct or dataset');
end

end