function menuH = createMenuItem(fig, parent, labelname, tag, varargin)
%CREATEMENUITEM Creates a menu item and connects its behavior
%   to the behavior of other ui components with specified tag name.
%
%   The returned menu item has 'tag_menu' as tag. The other specified UI
%   components can be toggle button or push button on toolbar.

% Copyright 2009 The MathWorks, Inc.


%== Create a menu item with parent.
if isempty(parent)
    menuH = uimenu(fig,'Label',labelname,'Tag', tag, varargin{:});
    return;
end
%== Create a menu item with tag_menu as name
menu_tag = [tag, '_menu'];
menuH = uimenu(parent,'Label',labelname,'Tag', menu_tag, varargin{:});

comp = findobj(fig,'tag',tag);
if ~isempty(comp)
    if strcmpi(get(comp,'type'),'uipushtool') % for toolbar pushbuttons
        set(menuH,'Callback',get(comp,'ClickedCallback'));
    elseif strcmpi(get(comp,'type'),'uitoggletool') % for toolbar toggles
        set(menuH,'Callback',get(comp,'ClickedCallback'));
    elseif strcmpi(get(comp,'type'),'uicontrol') %  for other push buttons
        set(menuH,'Callback',get(comp,'Callback'));
    end
end
end % end of function
