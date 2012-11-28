function x = peval(F,ds,J,nSets)
%function x = peval(F,ds,J,nSets)
% Evaluate x(:,i) = F{i}(ds) using parfor where F returns a matrix of double with
% the same number of rows as ds
%
% It is assumed that the calculation can be partitioned using a group
% assignment given by J. The vector J is a list of integer with the same
% length as size(ds,1)

if iscell(F),
   nCols = length(F);
else
   nCols = 1;
   F = {F};
end

n = size(ds,1);
assert(length(J)==n,'Expecting J to be a vector assigning an integer to each row');
assert(matlabpool('size')>0,'Expecting a pool of workers - use matlabpool to create one');
xI = cell(nSets);
I = equalPartition(nSets,J);
tI = cell(nSets,1);
for i=1:nSets,
    tI{i} = ds(I==i,:);
end
for colNo=1:nCols
    Fcol = F{colNo};
    parfor i = 1:nSets
        xI{i,colNo} = feval(Fcol,tI{i});
    end
end
x = nan(n,nCols);
for i=1:nSets,
    x(I==i,colNo) = xI{i,colNo};
end
end