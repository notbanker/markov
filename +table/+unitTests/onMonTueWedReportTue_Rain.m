function [dowName1,rain] = onMonTueWedReportTue_Rain(dowName,rain)
% An example of a user choice function
rain1 = rain;
switch dowName
    case {'M','T','W'},
        dowName1 = 'T';
    otherwise
        dowName1 = dowName;
end
end
