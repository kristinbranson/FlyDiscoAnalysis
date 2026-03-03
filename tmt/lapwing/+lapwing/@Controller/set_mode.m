function set_mode(self, desired_new_mode)
    % need to remember the old mode so that we can
    % uncheck that menu item
    old_mode = self.model.mode ;
    
    % set the chosen mode
    self.model.mode = desired_new_mode ;
    
    % Get the actual new mode
    new_mode = self.model.mode ;
    
    % unset the old mode button
    old_button_h=...
        findobj(self.figure_h,'Tag',sprintf('%s_button_h',old_mode));
    set(old_button_h,'Value',0);
    
    % set the new mode button
    new_button_h=...
        findobj(self.figure_h,'Tag',sprintf('%s_button_h',new_mode));
    set(new_button_h,'value',1);
end
