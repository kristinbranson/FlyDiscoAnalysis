function [allData,isNumericCol] = bioReadMixedData(fid,isNumericCol,delimiter,ignoreCols)
% BIOREADMIXEDDATA reads data from a table that contains mixed numeric and
% text columns. 
%
%   [ALLDATA, ISNUMERICCOL] = BIOREADMIXEDDATA(FID) reads tab delimited
%   data from columns of file handle FID and returns a cell array
%   containing the contents of the file and a logical array indicating
%   which columns contain only numeric data.
%
%   BIOREADMIXEDDATA(FID, ISNUMERICCOL) allows you to pass in a logical
%   array indicating the columns that may contain numeric data. You would
%   typically use this to indicate columns that are known to contain text
%   data.
%
%   BIOREADMIXEDDATA(FID, ISNUMERICCOL, DELIMITER) allows you to specify a
%   delimiter. Default is '\t'.
%
%   BIOREADMIXEDDATA(FID, ISNUMERICCOL, DELIMITER, IGNORECOLS) allows you
%   to specify a logical array of columns that should be ignored.
%
%   Examples:
%
%           [data, isnumericCols] = bioReadMixedData(fid,true(10,1);
%
%   See also TEXTSCAN.


% Copyright 2008 The MathWorks, Inc.


% Assume we are in a tab delimited w
if nargin < 3
    delimiter = '\t';
end

% if no information is given then guess everything is numeric
if nargin < 2
    pos = ftell(fid);
    line = fgetl(fid);
    tmpCols = strread(line,'%s','delimiter',delimiter);
    isNumericCol = true(1,numel(tmpCols));
    fseek(fid,pos,-1);
end

% Assume we want everything
if nargin < 4
    ignoreCols = false(size(isNumericCol));
end

% Create the format string
numCols = numel(isNumericCol);
formatStr = repmat('% s',1,numCols);
formatStr(3*find(ignoreCols)-1) = '*';
formatStr(3*find(isNumericCol)) = 'f';
% Keep columns that we want
colsToReadNdx = find(~ignoreCols);
% try to read everything in one go
loopNum = numCols;
fpos = ftell(fid);
while ~feof(fid) && loopNum > 0
    % note that we strip out spaces from format string
    allData = textscan(fid,strrep(formatStr,' ',''),'delimiter',delimiter);
    if ~feof(fid)
        % if textscan can't read the file, figure out which column
        % failed to be converted and then try again with this column as
        % a string
        badCol = find(diff(cellfun(@length,allData)))+1;
        isNumericCol(badCol) = false;
        formatStr(3*colsToReadNdx(badCol)) = 's';
        fseek(fid, fpos,-1);
        loopNum = loopNum -1; % should be redundant but let's be safe
    end
end
