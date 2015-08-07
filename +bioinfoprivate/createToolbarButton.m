function htool = createToolbarButton(toolbar, type, iconfile, tooltip, tag, varargin)
%CREATETOOLBARBUTTON Create toolbar buttons
% 
%  HT = CREATETOOLBARBUTTON(HTOOLBAR, TYPE, ICONFILE, TOOLTIP, TAG) returns
%  an UI button HT for toolbar HTOOLBAR. TYPE specifies the type of the
%  button. TYPE = 1 - uipushtool, 2-uitoggletool. ICONFILE specifies the
%  location of the icon image file.

% Copyright 2009 The MathWorks, Inc.


icon_cdata = bioinfoprivate.readIcon(iconfile);
props = {toolbar,'cdata', icon_cdata, 'Tooltip', tooltip, 'Tag', tag, varargin{:}}; 
switch type
    case 1
        htool=uipushtool(props{:});
    case 2
        htool=uitoggletool(props{:});
end
end % end of function
%--------------
