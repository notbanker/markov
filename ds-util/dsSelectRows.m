function dsNew = dsSelectRows(ds,I)
%function dsNew = dsSelectRows(ds,I)

if isa(ds,'dataset'),
    dsNew = dataset;
    vNames = ds.Properties.VarNames;
elseif isa(ds,'struct'),
    dsNew = struct;
    vNames = fieldnames(ds);
else
    error('Expecting struct or dataset');    
end

for k=1:length(vNames),
    vName = vNames{k};
    dsVName = ds.(vName);
    dsNew.(vName) = dsVName(I);
end


end