function [G,dsMin,In,New,Removed] = dependencyDag(function_names,ds)
% G = TABLE.PAR.DEPENDENCYDAG(FUNCTION_NAMES,DS) returns an adacency matrix
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

nFunc = length(function_names);
In = cell(nFunc,1);
New = cell(nFunc,1);
Removed = cell(nFunc,1);
dsMin = ds;
for k=1:nFunc,
    fn = function_names{k};
    [In{k},New{k},Removed{k}] = determineInOutRemoved(fn,ds);
    for m=1:length(New{k}),
        vn = New{k}{m};
        if ismember(vn,table.fieldnames(dsMin)),
            dsMin.(New{k}{m}) = [];
        end
    end
end
G = zeros(nFunc);
for i=1:nFunc,
    for j=1:nFunc
        if j~=i,
            G(i,j) = any(ismember(New{j},In{i})) | any(ismember(Removed{i},New{j})) | any(ismember(Removed{i},In{j}));
        end
    end
end
end

function [inFields,outFields,removedFields] = determineInOutRemoved(f,ds)
  missing = firstMissingVar(f,ds);
  while ~isempty(missing),
      disp(['Inspecting ',f,'. Adding the field ',missing,' filled with 666666 so topology can be determined']);
      ds.(missing) = 6666666*ones(size(ds,1),1);
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
    ds_without_fn = ds;
    fn = allArgs{m};
    ds_without_fn.(fn)=[];
    try
        throw_me_away = feval(f,ds_without_fn);
        dsMinimal.(fn) = [];
    catch ME
        % Leave the field in
    end
end
end


function fn = firstMissingVar(f,ds)
% Find the first missing variable name that is required to apply f to ds
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
    fn = token{1}{1};
 end
     

end

