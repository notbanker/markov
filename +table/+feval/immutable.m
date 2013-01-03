function [dsOut,ds,info] = immutable(function_names,ds,I,function_signatures)
% DSOUT = TABLE.FEVAL.IMMUTABLE(FUNCTION_NAMES,DS,I) takes a list of
% functions (or operators) that modify datasets or structs but treat fields as immutable
% in the sense that they do not modify existing fields entries but may delete
% the field outright. TABLE.FEVAL.IMMUTABLE determines a minimal dataset or
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
% [DSOUT,DS] = TABLE.PAR.FEVAL(FUNCTION_NAMES,DS,I) also returns
% the dataset or struct which is distributed to workers before an operator is applied. 
% This is the same as ds but with some fields removed. 
%
% See also parfor_ds_sequential

%TODO:Treat cycles in graphs gracefully or warn
%TODO:More graceful exception handling 

if ischar(function_names),
    function_names = {function_names}; 
end

% Use a small sample to determine topological ordering and minimal input arguments
dsTest = table.select.rows(ds,I==1);
if table.size(dsTest,1)>5,
   dsTest = table.select.rows(dsTest,1:5); 
end
if nargin>=4,
    [G,dsMin,info.function_signatures] = dependencyGraph(function_names,dsTest,function_signatures);
else
    [G,dsMin,info.function_signatures] = dependencyGraph(function_names,dsTest);
end
info.seq = toposort(G');
info.orderedFunctions = function_names(info.seq)';
info.minimalFields = table.fieldnames(dsMin);

% Create a smaller dataset or struct to send to workers
fns = table.fieldnames(ds);
fns_ = table.fieldnames(dsMin);
for k=1:length(fns),
   if ~ismember(fns{k},fns_),
      ds = table.rmfield(ds,fns{k}); 
   end
end

% Warn if we need to add more variables
fns = table.fieldnames(ds);
for k=1:length(fns),
   if all(666666==dsMin.(fns{k})),
      error(['Please add the field ',fns{k},' to the dataset or struct']);    
   end
end

% Apply parfor_ds_sequential in sensible order with reduced overhead
isImmutable = true;
tic,
dsOut = table.feval.sequential(function_names(info.seq),ds,I,isImmutable);
info.cpu = toc;
end

function [G,dsMin,function_signatures] = dependencyGraph(function_names,ds,function_signatures)
% G = TABLE.PAR.DEPENDENCYGRAPH(FUNCTION_NAMES,DS) returns an adacency matrix
% G where G(i,j)==1 if the j'th function needs to be computed before the i'th
% function, either because the i'th function requires a field that is an output
% of the j'th function; or because the j'th function removes a field that is an 
% output of the i'th function; or because the i'th function removes a field
% that is an input to the j'th function. 
%
% DS should be a maximal 'test' dataset that contains
% the union of all fields that are required by functions in the list
%
% Though all functions modify datasets it is presumed that they treat
% fields as immutable. Thus their action is the addition of new fields. 

%% Determine signatures if we need to
nFunc = length(function_names);
n = table.size(ds,1);
if nargin>=3,
    assert(isa(function_signatures,'struct'),'Expecting function_signatures to be a struct');
    assert(all(table.isfield(function_signatures,{'In','New','Removed'})),'Expecting fields In, New and Removed');
elseif nargin<3,
    function_signatures = struct;
    function_signatures.In = cell(nFunc,1);
    function_signatures.New = cell(nFunc,1);
    function_signatures.Removed = cell(nFunc,1);
    for k=1:nFunc,
        fn = function_names{k};
        [function_signatures.In{k},function_signatures.New{k},function_signatures.Removed{k}] = determineInOutRemoved(fn,ds);
    end
end

%% Determine minimal table to pass
dsMin = struct;
for k=1:nFunc
    for m=1:length(function_signatures.In{k}),
       fn = function_signatures.In{k}{m};
       if ~table.isfield(dsMin,fn),
           if table.isfield(ds,fn),
               dsMin.(fn) = ds.(fn);
           else
               dsMin.(fn) = 666666*ones(n,1);
           end
       end
    end
end
for k=1:nFunc,
    for m=1:length(function_signatures.New{k}),
        fn = function_signatures.New{k}{m};
        if table.isfield(dsMin,fn),
            dsMin = table.rmfield(dsMin,function_signatures.New{k}{m});
        end
    end
end

%% Construct dependency graph
G = zeros(nFunc);
for i=1:nFunc,
    for j=1:nFunc
        if j~=i,
            G(i,j) = any(ismember(function_signatures.New{j},function_signatures.In{i}))...
                  | any(ismember(function_signatures.Removed{i},function_signatures.New{j})) | ...
                  any(ismember(function_signatures.Removed{i},function_signatures.In{j}));
        end
    end
end
end

function [inFields,outFields,removedFields] = determineInOutRemoved(f,ds)
  missing = firstMissingVar(f,ds);
  while ~isempty(missing),
      %disp(['While inspecting ',f,' added the field ',missing,' filled with 666666 so topology can be determined']);
      ds.(missing) = 666666*ones(table.size(ds,1),1);
      missing = firstMissingVar(f,ds);
  end
  dsMinimal = removeUnusedInputFields(f,ds);
  dsOut = feval(f,dsMinimal);
  inFields = table.fieldnames(dsMinimal);
  outFields = setdiff(table.fieldnames(dsOut),table.fieldnames(dsMinimal));
  removedFields = setdiff(table.fieldnames(dsMinimal),table.fieldnames(dsOut));
end


function dsMinimal = removeUnusedInputFields(f,ds)
dsMinimal = ds;
allArgs = table.fieldnames(ds);
nArgs = length(allArgs);
for m=1:nArgs,
    fn = allArgs{m};
    ds_without_fn = table.rmfield(ds,fn);
    try
        throw_me_away = feval(f,ds_without_fn);
        dsMinimal = table.rmfield(dsMinimal,fn);
    catch ME
        % Leave the field in
    end
end
end


function fn = firstMissingVar(f,ds)
% Find the first missing variable name that is required to apply f to ds
% FIXME: Very brittle. Assumes that all operators fail only because they
% call a missing field. 
 fn = '';
 try
    feval(f,ds);
 catch ME
    msg = ME.message;
    if isa(ds,'dataset'),
        pat = 'Unrecognized variable name ''([\w\d]+)';
    else
        pat = 'Reference to non-existent field ''([\w\d]+)';
    end
    [~,~,~,~,token] = regexp(msg, pat);
    if isempty(token),
        warning(['firstMissingVar failed because operator ',f,' returned an error other than a missing field: ',msg]);
        error('Not sure what to do');
    else
        fn = token{1}{1};
    end
 end
     

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


