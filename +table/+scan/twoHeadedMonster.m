function ds = twoHeadedMonster(filename,varargin)
% DS = TABLE.SCAN.TWOHEADEDMONSTER(FILENAME) reads potentially huge 
% dataset from a delimited file with two header lines. 
%
% TABLE.SCAN.TWOHEADEDMONSTER uses a lousy, linux depenent hack to break up
% huge files and pass them to TABLE.SCAN.TWOHEADED
%
% The first header line specifies variable names. 
% The second header line specifies data types
% 
% Data types are belonging to one of the following:
%  'double'    Floating point
%  'string'    Strings of possibly unequal length
%  'char'      Strings of equal length
%  'int'       Integers
%  'datetime'  Char that can be converted to datetime
%
% If a column is of type datetime then the first data entry, on the 3rd
% line of the file, is used to infer the date format used. 


if ~isunix,
   error('Only written for linux. If your files are small enough try table.scan.twoHeaded instead.'); 
end

if isempty(strfind(filename,filesep)),
   filename = which(filename);
end

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

nLines = countLinesInAFileEfficiently(filename);

if (nLines<varargin{3}+500)
   % If the file is small enough then read it in directly
   ds = table.scan.twoHeaded(filename,varargin{:});
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
   dsRest = table.scan.twoHeadedMonster(rest,varargin{:});
   dsStart = table.scan.twoHeaded(start,varargin{:});
   system(['rm ',rest]);
   system(['rm ',start]);
   ds = [dsStart;dsRest];
end

end


function n = countLinesInAFileEfficiently(filename)
fh = fopen(filename, 'r');
chunksize = 1e6; % read chuncks of 1MB at a time
n = 0;
while ~feof(fh)
    ch = fread(fh, chunksize, '*uchar');
    if isempty(ch)
        break
    end
    n = n + sum(ch == sprintf('\n'));
end
fclose(fh);
end