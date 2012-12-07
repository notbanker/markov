function ds = twoscan(filename,varargin)
% DS = TWOSCAN(FILENAME) reads a dataset from a
%
% Read dataset from delimited file with two line header. The first
% line is the variable name. The second is the type, belonging to 
% {double, datetime, string, char, int}

filename = GetFullPath(filename);

%maxlines
if nargin<4
   varargin{3} = 1000000; 
end

% default is to convert dates to matlab datenum
if nargin<3,
   varargin{2} = true; 
end

% default delimiter is .csv
if nargin<2,
    varargin{1} = ',';
end

nLines = countLinesInAFile(filename);
if (nLines<varargin{3}+500)
   % If the file is small enough then read it in directly
   ds = twoscan1(filename,varargin{:});
else
   % A bit clunky: break into two problems and use the file system for tail recursion
   tmpHeader = [filename,'-tmp-header'];
   start  = [filename,'1'];
   tmpRest   = [filename,'-tmp-rest'];
   start = [filename,'1'];
   rest = [filename,'2'];
   system(['head -n 3 ',filename,' > ',tmpHeader]);
   system(['head -n ',num2str(varargin{3}),' ',filename,' > ',start]);
   system(['tail -n +',num2str(varargin{3}+1),' ',filename,' > ',tmpRest]);
   system(['cat ',tmpHeader,' ',tmpRest,' > ',rest]);
   system(['rm ',tmpHeader]);
   system(['rm ',tmpRest]);
   dsRest = twoscan(rest,varargin{:});
   dsStart = twoscan1(start,varargin{:});
   system(['rm ',rest]);
   system(['rm ',start]);
   ds = [dsStart;dsRest];
end

end

function ds = twoscan1(filename,delimiter,convertToDatenum,dummy)

% textscan formats
formats.default = '%s';
formats.double = '%f';
formats.datetime = '%s';
formats.string = '%s';  
formats.char = '%s';
formats.int = '%d';

% headers & first example
try
    fid = fopen(filename,'rt');
    header1 = fgetl(fid);
    header2 = fgetl(fid);
    line3 = fgetl(fid);
    fclose(fid);
catch ME
    error(['File must contain at least three lines ',ME.message]);
end

% get variable names
splits1 = [0,find(header1 == delimiter), length(header1)+1];
fields = cell(length(splits1)-1,1);
for i = 1:length(splits1)-1
	fields{i} = strtrim(header1(splits1(i)+1:splits1(i+1)-1));
end  

% get type names
splits2 = [0,find(header2 == delimiter), length(header2)+1];
assert(length(splits2)==length(splits1),['There are ',num2str(length(splits1)),...
        ' variable names on first header line but ',num2str(length(splits2)),' type names listed on the second header line.']);
types = cell(length(splits2)-1,1);
isDate = false(length(splits2),1);
isString = false(length(splits2),1);
for i = 1:length(splits2)-1
	types{i} = lower(strtrim(header2(splits2(i)+1:splits2(i+1)-1)));
    isDate(i) = strcmpi(types{i},'datetime');
    isString(i) = strcmpi(types{i},'string');
end  

% get first example
splits3 = [0,find(line3 == delimiter), length(line3)+1];
assert(length(splits3)==length(splits1),['There are ',num2str(length(splits1)),...
        ' variable names on first header line but ',num2str(length(splits3)),' entries on the third row of the file.']);
examples = cell(length(splits3)-1,1);
for i = 1:length(splits3)-1
	examples{i} = strtrim(line3(splits3(i)+1:splits3(i+1)-1));
end  

% interpret as textscan formats
textFormats = cell(length(fields), 1);
for i = 1:length(fields)
    if isfield(formats,types{i}),
        textFormats{i} = formats.(types{i});
    else
        textFormats{i} = formats.default;
    end
end

% textscan the file
formatString = cell2mat(textFormats');
fid = fopen(filename);
d = textscan(fid, formatString, 'headerlines', 2, 'delimiter', delimiter, 'TreatAsEmpty', 'null');
fclose(fid);

% convert to dataset
nRows = length(d{1});
nCols = length(d);
if isequal(nRows, 0),
    ds = [];
else
    data = cell(nRows, nCols);
    for i = 1:nCols
        if isequal(class(d{i}),'cell'),
            data(:,i) = d{i};
        else
            data(:,i) = num2cell(d{i});
        end
    end
end

% convert to dataset and store text and data formats in UserData
if ~isequal(nRows,0),
    ds = dataset({data, fields{:}});
    for i=1:length(fields),
        if ~isString(i),
            fn = fields{i};
            ds.(fn) = cell2mat(ds.(fn));
        end
    end
    
    % Determine datetime formats and convert
    ds.Properties.UserData.textFormats = textFormats;
    dateFormats = cell(nCols,1);
    for colNo =1:nCols,
        if isDate(colNo),
            eg = examples{colNo};
            dateFormats{colNo} = getDateFormat(eg);
            if convertToDatenum,
                ds.(fields{colNo}) = datenum(ds.(fields{colNo}),dateFormats{colNo});
                types{colNo} = 'datenum';
            end
        end
    end
    ds.Properties.UserData.varNames = ds.Properties.VarNames';
    ds.Properties.UserData.varTypes = types;
    ds.Properties.UserData.varDateFormats = dateFormats;
end

end

function [format] = getDateFormat(date_string)
% Ripped from datevec and modified to include our formats

formats = {'yyyy-mm-dd HH:MM:SS',...
    'yyyy-mm-dd',...
    'dd-mmm-yyyy HH:MM:SS',...
    'mm/dd/yy',...
    'mm/dd',...
    'HH:MM:SS',...
    'HH:MM:SS PM',...
    'HH:MM',...
    'HH:MM PM',...
    'mm/dd/yyyy',...
    'dd-mmm-yyyy HH:MM',...
    'dd-mmm-yy'};
nFormats = length(formats);

i = 1;
foundIt = false;
while i<nFormats && ~foundIt,
    try
        D = datenum(date_string,formats{i});
        D1 = datenum(datestr(D,formats{i}));
        if D == D1,
            format = formats{i};
            foundIt = true;
        end
    catch
        i = i+1;
    end
end

end


