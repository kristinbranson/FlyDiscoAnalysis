function state = toggleState(hfig, hSrc)
%TOGGLESTATE Return and set a toggele button state and related menu item state.
%
% Useful for syncing toolbar toggle button with a menu item. The menu item
% should be created by the bioinfoprivate.creatMenuItem function. 

% Copyright 2009 The MathWorks, Inc.


%== Figure out hSrc type.
type = get(hSrc,'type');
if strcmpi(type,'uimenu')
    tag = strrep(get(hSrc,'tag'),'_menu','');
    hSrcPeer = findall(hfig,'tag', tag);

    if strcmp(get(hSrcPeer,'State'),'on')
        set(hSrcPeer,'State','off');
    else
        set(hSrcPeer,'State','on');
    end
    state = get(hSrcPeer, 'State');
    set(hSrc, 'Checked', state);
elseif strcmpi(type, 'uitoggletool')
    menu_tag = [get(hSrc,'tag'),'_menu'];
    hSrcPeer = findall(hfig, 'tag', menu_tag);
    state = get(hSrc,'State');

    if ~isempty(hSrcPeer)
        set(hSrcPeer,'Checked', state);
    end
end
end % end of function
