function mouse_button_down_in_image(self)
    
    % get the current mode
    mode = self.model.mode ;
    
    % switch on the current roi mode appropriately
    switch(mode)
        case 'zoom'
            selection_type = get(self.figure_h,'SelectionType') ;
            switch selection_type
                case 'extend'
                    self.zoom_out();
                case 'alternate'
                    self.zoom_out();
                case 'open'
                    self.zoom_out();
                otherwise
                    self.draw_zoom_rect('start');
            end
    end
    
end
