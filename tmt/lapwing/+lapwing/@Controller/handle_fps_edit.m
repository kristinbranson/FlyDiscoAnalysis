function handle_fps_edit(self)

new_fps_string=get(self.FPS_edit_h,'String');
new_fps=str2double(new_fps_string);
if isfinite(new_fps) && (new_fps>0)
  self.model.fps=new_fps;
end
self.update_fps();

end
