function quit(self)
    self.stop_playing();
    pause(0.01);
    self.model.close_file();
    self.close();
end
