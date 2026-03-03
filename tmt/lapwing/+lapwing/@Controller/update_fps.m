function update_fps(self)
    FPS_edit_string = sprintf('%g', self.model.fps) ;
    set(self.FPS_edit_h, 'String', FPS_edit_string) ;
end
