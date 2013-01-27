function thisDirectory = markov
% Returns top level directory where this file resides
    if (~isdeployed)
        w = which('branch');
        thisDirectory = w(1:end-length('branch.m')-1);
    else
        thisDirectory = pwd;
    end
end