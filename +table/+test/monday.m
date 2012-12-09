function [dowName1,rain1] = monday(dowName,rain)
switch lower(dowName)
    case 's'
        dowName1 = 'F';
        rain1 = rain+1;
    otherwise
        dowName1 = 'T';
        rain1 = rain-1;
end
end