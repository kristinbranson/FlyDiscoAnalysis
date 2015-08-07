function rgbVal = colorSpecLookUp(colName)
%COLORSPECLOOKUP Return the RGB value of predefined colors.
%
%   RGBVAL = COLORSPECLOOKUP(COLNAME) returns the RGB value of eight
%   predefined colors: yellow(y), magenta(m), cyan(c), red(r), green(g),
%   blue(b), white(w), black(k). See ColorSpec documentation for more
%   details. COLNAME must be the predfine color names.

% Copyright 2009 The MathWorks, Inc.


rgbVal = [ 0 0 0];
%==Define color names according to ColorSpec
shortNames = {'y', 'm', 'c', 'r', 'g', 'b', 'w', 'k'};
longNames = {'yellow', 'magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black'};
RGBS = [ 1 1 0;
         1 0 1;
         0 1 1;
         1 0 0;
         0 1 0;
         0 0 1;
         1 1 1;
         0 0 0];
%== Lookup
k = find(strncmpi(colName, [shortNames,longNames], numel(colName)), 1);
if ~isempty(k)
    if k > 8
       k = k - 8;
   end
   rgbVal = RGBS(k, :);
else
   bioinfoprivate.bioerror(mfilename, 'WrongColorName',...
        'The input color name did not match the predefined color names.')
end
end
