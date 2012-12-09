function [dowName1] = onMonTueWedReportTue(dowName)
% An example of a user choice function
switch dowName
    case {'M','T','W'},
        dowName1 = 'T';
    otherwise
        dowName1 = dowName;
end
end