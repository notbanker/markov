function dsNew = cols(ds,field_names)
%function dsNew = cols(ds,field_names)

old_names = table.fieldnames(ds);
if isa(ds,'dataset'),
    dsNew = dataset;
elseif isa(ds,'struct'),
    dsNew = struct;
end

for k=1:length(field_names),
   fn = field_names{k};
   if ismember(fn,old_names),
      dsNew.(fn) = ds.(fn); 
   else
       warning(['Attempt to select ds.(',fn,') which is not a field']);
   end
end


end