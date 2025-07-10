function set_cmap_name(self, new_cmap_name)
    % need to remember the old cmap_name so that we can
    % uncheck that menu item
    old_cmap_name = self.model.cmap_name ;

    % set the chosen cmap_name
    self.model.cmap_name = new_cmap_name ;
    
    % Get the new colormap
    cmap = self.model.cmap ;

    % uncheck the old menu item
    menu_h=eval(sprintf('self.%s_menu_h',old_cmap_name));
    set(menu_h,'Checked','off');

    % check the new menu item
    menu_h=eval(sprintf('self.%s_menu_h',new_cmap_name));
    set(menu_h,'Checked','on');

    % set the colormap
    set(self.figure_h,'colormap',cmap);
end
