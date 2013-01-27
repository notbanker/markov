function thisDirectory = markov
% Returns top level directory where this file resides
    if (~isdeployed)
        w = which('markov');
        thisDirectory = w(1:end-length('markov.m')-1);
    else
        thisDirectory = pwd;
    end
end