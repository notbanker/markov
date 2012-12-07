function [dsOut,ds] = parfor_ds_immutable(function_names,ds,I)
% DSOUT = PARFOR_DS_IMMUTABLE(FUNCTION_NAMES,DS,I) takes a list of
% functions (or operators) that modify datasets or structs but treat fields as immutable
% in the sense that they do not modify existing fields entries but may delete
% the field outright. PARFOR_DS_IMMUTABLE determines a minimal dataset or
% struct dsIn and then sequentially applies the operators, not necessarily in the
% supplied order, but rather in a topological ordering as determined by the
% needs of the operators (some operators require input fields that are an
% output of a previous operation, and this defines the topological ordering). 
% 
% Operations are parallelized on row groups. It is assumed that the calculation
% can be partition using a row group assignment given by I. The vector I is a
% list of integers with the same length as dsSize(ds,1). 
%
% An easy way to create the partition is calling I = equalPartition(nGroups,J1,J2,...) 
% This will divide the rows into roughly equal sized groups in such a way
% that all rows in a group take the same values for J1, J2. In particular
% the J1, J2, ... might be fields of ds taking discrete values. 
%
% [DSOUT,DS] = PARFOR_DS_IMMUTABLE(FUNCTION_NAMES,DS,I) also returns
% the dataset or struct which is distributed to workers before an operator is applied. 
% This is the same as ds but with some fields removed. 
%
% See also parfor_ds_sequential

if ischar(function_names),
    function_names = {function_names}; 
end

% Use a small sample to determine topological ordering and minimal input arguments
dsTest = dsSelectRows(ds,I==1);
if dsSize(dsTest,1)>50,
   dsTest = dsSelectRows(dsTest,1:50); 
end
[G,dsMin] = ds_dependency_dag(function_names,dsTest);
seq = toposort(G');

% Create a smaller dataset or struct to send to workers
fns = dsNames(ds);
fns_ = dsNames(dsMin);
for k=1:length(fns),
   if ~ismember(fns{k},fns_),
      ds.(fns{k}) = []; 
   end
end

% Warn if we need to add more variables
fns = dsNames(ds);
for k=1:length(fns),
   if all(666666==dsMin.(fns{k})),
      error(['Please add the field ',fns{k},' to the dataset or struct']);    
   end
end

% Apply parfor_ds_sequential in sensible order with reduced overhead
immutable = true;
dsOut = parfor_ds_sequential(function_names(seq),ds,I,immutable);

end


function [seq] = toposort(adj)
% TOPOSORT		A Topological ordering of nodes in a directed graph
%
%  [SEQ] = TOPOSORT(ADJ)
%
% Inputs :
%    ADJ : Adjacency Matrix.
%	   ADJ(i,j)==1 ==> there exists a directed edge
%	   from i to j
%
% Outputs :
%    SEQ : A topological ordered sequence of nodes.
%          empty matrix if graph contains cycles.
%
% Usage Example :
%		N=5;
%		[l,u] = lu(rand(N));
%		adj = ~diag(ones(1,N)) & u>0.5;
%		seq = toposort(adj);
%
%
% Note     :
% See also

% Uses :

% Change History :
% Date		Time		Prog	Note
% 18-May-1998	 4:44 PM	ATC	Created under MATLAB 5.1.0.421

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl

N = size(adj);
indeg = sum(adj,1);
outdeg = sum(adj,2);
seq = [];

for i=1:N,
    % Find nodes with indegree 0
    idx = find(indeg==0);
    % If can't find than graph contains a cycle
    if isempty(idx),
        seq = [];
        break;
    end;
    % Remove the node with the max number of connections
    [~, idx2] = max(outdeg(idx));
    indx = idx(idx2);
    seq = [seq, indx];
    indeg(indx)=-1;
    idx = find(adj(indx,:));
    indeg(idx) = indeg(idx)-1;
end
end


