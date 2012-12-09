function names = fieldnames(ds)

if isa(ds,'dataset'),
    names = ds.Properties.VarNames;
elseif isa(ds,'struct'),
    names = fieldnames(ds);
end
    
end