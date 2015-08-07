
function property = biofindprop(obj, propname)
%BIOFINDPROP  Find a property. 
%
%   P = BIOFINDPROP(OBJ, PROPNAME) for finds and returns the property
%   object associated with property name PROPNAME of a handle object OBJ
%   for HG1 or HG2.
 
%   Copyright 2008 The MathWorks, Inc.

 
% HG1/HG2 Safe way find a property
if feature( 'HGUsingMATLABClasses' )
    property = findprop(obj, propname);
else
    % Convert double to handle is allowed because this branch is under HG1
    ws = warning( 'off', 'MATLAB:hg:DoubleToHandleConversion' );
    obj = handle( obj );
    warning( ws );
    % Make listener
    property = findprop(obj, propname);
end
 
end % biofindprop
