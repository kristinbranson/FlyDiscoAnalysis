function handle_key_press(self,event)

% handle the key
if strcmp(event.Character,',') || strcmp(event.Character,'<')
  self.change_z_slice_rel(-1);
elseif strcmp(event.Character,'.') || strcmp(event.Character,'>')
  self.change_z_slice_rel(+1);
elseif strcmp(event.Character,'p')    
  self.play(+1);
elseif strcmp(event.Key,'delete') || strcmp(event.Key,'backspace')
  self.delete_selected_roi();  
end

% % update the shift state, if needed
% if strcmp(event.Key,'shift')
%   %fprintf('shift pressed\n');
%   self.shift_depressed=true;
% end

end
