function ds = sequential(function_names,ds,I,immutable)
% DS = TABLE.FEVAL.SEQUENTIAL(FUNCTION_NAMES,DS,I) sequentially applies a list
% of functions in cell array of char FUNCTION_NAMES to a dataset ds, but parallelizes
% the operations into groups of rows.  
% 
% It is assumed that each function in the list of FUNCTION_NAMES takes a
% dataset or tabular struct as argument and returns another dataset. Each function in
% FUNCTION_NAMES finds interpretation as a bulk update acting on multiple
% records (i.e. rows) and possibly creating new fields. 
%
% It is assumed that the calculation can be partitioned using a row group
% assignment given by I. The vector I is a list of integers with the same
% length as table.size(ds,1). 
%
% An easy way to create the partition is calling I = equalPartition(nGroups,J1,J2,...) 
% This will divide the rows into roughly equal sized groups in such a way
% that all rows in a group take the same values for J1, J2. In particular
% the J1, J2, ... might be fields of ds taking discrete values. 
%
% DS = TABLE.FEVAL.SEQUENTIAL(FUNCTION_NAMES,DS,I,true) instructs table.feval.sequential that
% fields are immutable. If this flag is set the only fields that will change in 
% ds are those that are either newly created or deleted outright by the
% cumulative action of the functions in the list. 
%
% See also table.feval.immutable

nGroups = max(I);

if nargin<4,
   immutable = false;
   warning('You might increase performance by specifing that fields are immutable');
end

if ~isa(function_names,'cell'),
    function_names = {function_names}; 
end

n = table.size(ds,1);
assert(length(I)==n,'Expecting I to be a vector assigning an integer to each row');
assert(matlabpool('size')>0,'Expecting a pool of workers - use matlabpool to create one');
dsOut = cell(nGroups,1);
dsIn = cell(nGroups,1);
for i=1:nGroups,
    dsIn{i} = table.selectRows(ds,I==i);
end

nFunc = length(function_names);
parfor i = 1:nGroups
    ds_ = dsIn{i};
    for k=1:nFunc,
        ds_ = feval(function_names{k},ds_);
    end
    dsOut{i} = ds_;
end

newFields = setdiff(table.fieldnames(dsOut{1}),table.fieldnames(dsIn{1}));
if immutable,
    potentiallyChangedFields = newFields;
else
    potentiallyChangedFields = table.fieldnames(dsOut{1}); % all of them
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

removedFields = setdiff(table.fieldnames(dsIn{1}),table.fieldnames(dsOut{1}));
nRemoved = length(removedFields);
for k=1:nRemoved,
   fn = removedFields{k};
   ds.(fn) = [];
end

end