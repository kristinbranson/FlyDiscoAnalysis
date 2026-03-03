function change_z_slice_abs(self, desired_new_z_index)
    % Change the current z_slice to the given z_slice index.
    
    old_z_index = self.model.z_index ;
    self.model.z_index = desired_new_z_index ;
    new_z_index = self.model.z_index ;
    if new_z_index ~= old_z_index ,
        set(self.z_index_edit_h, 'String', sprintf('%d',new_z_index)) ;
        set(self.image_h, 'CData', self.model.indexed_z_slice) ;
    end
end
