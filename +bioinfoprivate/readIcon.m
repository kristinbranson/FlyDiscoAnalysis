function cdata = readIcon(filename)
%READICON Read an icon image file and convert it to CData.
%
% Useful for reading toolbar icon image files in true color.

% Copyright 2009 The MathWorks, Inc.


[p,f,ext] = fileparts(filename); 
% if this is a mat-file, look for the variable cdata 
if isequal(lower(ext),'.mat')
    data = load(filename,'cdata');
    cdata = data.cdata;
    return;
end

[cdata,map] = imread(filename);
if isempty(cdata)
    return;
end

if isempty(map)
    % need to use doubles, nan's only work as doubles
    cdata = double(cdata);
    cdata = cdata/255;
else
    % Set all white (1,1,1) colors to be transparent (nan)
    ind = find(map(:,1)+map(:,2)+map(:,3)==3);
    map(ind) = nan; %#ok
    cdata = ind2rgb(cdata,map);
end
end
