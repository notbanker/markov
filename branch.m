function thisDirectory = branch
% Put this at the top level of the modules directory (synomym for branch)
    if (~isdeployed)
        w = which('branch');
        thisDirectory = w(1:end-length('branch.m')-1);
    else
        thisDirectory = pwd;
    end
end