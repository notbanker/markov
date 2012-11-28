function ds = parfor_dataset_sequential(function_names,ds,I,immutable)
% DS = PARFOR_DATASET_SEQUENTIAL(FUNCTION_NAMES,DS,I) sequentially applies a list
% of functions to a dataset ds, but parallelizes the operations into groups
% of rows.  
% 
% It is assumed that each function in the list of FUNCTION_NAMES takes a
% dataset as argument and returns another dataset. 
%
% It is assumed that the calculation can be partition using a row group
% assignment given by I. The vector I is a list of integers with the same
% length as size(ds,1). 
%
% An easy way to create the partition is calling I = equalPartition(nGroups,J1,J2,...) 
% This will divide the rows into roughly equal sized groups in such a way
% that all rows in a group take the same values for J1, J2. In particular
% the J1, J2, ... might be fields of ds taking discrete values. 
%
% DS = PARFOR_DATASET_SEQUENTIAL(FUNCTION_NAMES,DS,I,true) instructs parfor_dataset_sequential that
% fields are immutable. If this flag is set the only fields that will change in 
% ds are those that are either newly created or deleted outright by the
% cumulative action of the functions in the list. 
%
% See also parfor_dataset_immutable

nGroups = max(I);

if nargin<4,
   immutable = false;
   warning('You might increase performance by specifing that fields are immutable');
end

if ~isa(function_names,'cell'),
    function_names = {function_names}; 
end

n = size(ds,1);
assert(length(I)==n,'Expecting I to be a vector assigning an integer to each row');
assert(matlabpool('size')>0,'Expecting a pool of workers - use matlabpool to create one');
dsIn = cell(nGroups,1);
dsOut = cell(nGroups,1);
for i=1:nGroups,
    dsIn{i} = ds(I==i,:);
end
nFunc = length(function_names);
parfor i = 1:nGroups
    ds_ = dsIn{i};
    for k=1:nFunc,
        ds_ = feval(function_names{k},ds_);
    end
    dsOut{i} = ds_;
end

newFields = setdiff(dsOut{1}.Properties.VarNames,dsIn{1}.Properties.VarNames);
if immutable,
    potentiallyChangedFields = newFields;
else
    potentiallyChangedFields = dsOut{1}.Properties.VarNames; % all of them
end

nPotentiallyChanged = length(potentiallyChangedFields);
for k=1:nPotentiallyChanged,
    fn = potentiallyChangedFields{k};
    x = nan(n,size(dsOut{1}.(fn),2));
    for i=1:nGroups,
        x(I==i,:) = dsOut{i}.(fn);
    end
    ds.(fn) = x;
end

removedFields = setdiff(dsIn{1}.Properties.VarNames,dsOut{1}.Properties.VarNames);
nRemoved = length(removedFields);
for k=1:nRemoved,
   fn = removedFields{k};
   ds.(fn) = [];
end

end