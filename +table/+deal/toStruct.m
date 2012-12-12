function dsStruct = toStruct(dsDataset)
% Memory efficient conversion of dataset to struct
if isstruct(dsDataset),
    dsStruct = dsDataset;
elseif isa(dsDataset,'dataset')
    vn = table.fieldnames(dsDataset);
    dsStruct = struct;
    for k=1:length(vn),
        dsStruct.(vn{k}) = dsDataset.(vn{k});
        dsDataset.(vn{k}) = []; % conserve mem as we go
    end
end

end