function ds = peval_struct(F,ds,I)
%
% Adds new fields using parfor by sequentially applying one or more
% functions F{k}. 
%
% F{k} takes struct or dataset ds and returns the same with more fields
%
% It is assumed that the calculation can be partition using a group
% assignment given by I. The vector I is a list of integers with the same
% length as size(ds,1). Use I = equalPartition(nSets,J) to create an
% assignment using some non-continguous index J.


nSets = max(I);

wasStruct = isa(ds,'struct');
if wasStruct,
   ds = struct2dataset(ds); 
end

if isa(F,'function_handle'),
    F = {F}; 
end

n = size(ds,1);
assert(length(I)==n,'Expecting I to be a vector assigning an integer to each row');
assert(matlabpool('size')>0,'Expecting a pool of workers - use matlabpool to create one');
dsIn = cell(nSets,1);
dsOut = cell(nSets,1);
for i=1:nSets,
    dsIn{i} = ds(I==i,:);
end

parfor i = 1:nSets
    ds_ = dsIn{i};
    for k=1:nF
        ds_ = feval(F{k},ds_);
    end
    dsOut{i} = ds_;
end

newFields = setdiff(dsOut{1}.Properties.VarNames,dsIn{1}.Properties.VarNames);

nNew = length(newFields);
for k=1:nNew,
    fn = newFields{k};
    x = nan(n,size(dsOut{1}.(fn),2));
    for i=1:nSets,
        x(I==i,:) = dsOut{i}.(fn);
    end
    ds.(fn) = x;
end

if wasStruct,
   ds = dataset2struct(ds); 
end

end