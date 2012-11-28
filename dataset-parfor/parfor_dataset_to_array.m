function [x,start_ndx,end_ndx] = parfor_dataset_to_array(function_names,ds,I)
%function [x,start_ndx,end_ndx] = parfor_dataset_to_array(function_names,ds,I)
%
% X = PEVAL_TO_ARRAY(FUNCTION_NAMES,DS,I) concatenates columnwise the results of
% several functions, parallelizing operations by groups of rows. Each function takes
% the dataset DS as argument and returns a matrix with the same number of rows as DS. 
%
% [X,START_NDX,END_NDX] = PEVAL_TO_ARRAY(FUNCTION_NAMES,DS,I) also returns the starting and
% ending column indexes for the result of function_names{i} 
%
% It is assumed that the calculation can be partitioned using a group
% assignment given by I. The vector I is a list of integers with the same
% length as size(ds,1).
%
% An easy way to create the partition is calling I = equalPartition(nGroups,J1,J2,...) 
% This will divide the rows into roughly equal sized groups in such a way
% that all rows in a group take the same values for J1, J2. In particular, one might
% set J1, J2,... equal to one or more fields of DS provided that the fields take on a finite
% number of values only. 

if iscell(function_names) 
   nFunc = length(function_names);
else
   nFunc = 1;
   function_names = {function_names};
end

n = size(ds,1);
assert(length(I)==n,'Expecting I to be a vector assigning an integer to each row');
assert(matlabpool('size')>0,'Expecting a pool of workers - use matlabpool to create one');
xI = cell(nGroups);
tI = cell(nGroups,1);
for i=1:nGroups,
    tI{i} = ds(I==i,:);
end
for funcNo=1:nFunc
    Fcol = function_names{funcNo};
    parfor i = 1:nGroups
        xI{i,funcNo} = feval(Fcol,tI{i});
    end
end
start_ndx = nan(nFunc,1);
end_ndx = nan(nFunc,1);
offset = 0;
for funcNo=1:nFunc,
   len = size(xI{1,funcNo},2);
   start_ndx(funcNo) = offset+1;
   end_ndx(funcNo) = offset+len;
   offset = offset+len;
end

x = nan(n,offset);
for i=1:nGroups,
    x(I==i,start_ndx(funcNo):end_ndx(funcNo)) = xI{i,funcNo};
end
end