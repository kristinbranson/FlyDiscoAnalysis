function update(self)
    
    % Called when the model changes.  Updates the view to reflect the current
    % model.  This is the 'nuclear option' of update functions.  This gets all
    % aspects of the model that affect the view, and updates those parts of the
    % view appropriately.
    
    % get the local filename
    [~,base_name,ext]=fileparts(self.model.file_name);
    file_name_local=[base_name ext];
    
    % set the figure title
    if isempty(file_name_local)
        title_string='Lapwing';
    else
        title_string=sprintf('Lapwing - %s',file_name_local);
    end
    set(self.figure_h,'name',title_string);
    
    % % determine the colorbar bounds
    % [data_min,data_max]=self.model.pixel_data_type_min_max();
    % self.colorbar_min_string=sprintf('%d',data_min);
    % self.colorbar_max_string=sprintf('%d',data_max);
    % self.colorbar_min=str2double(self.colorbar_min_string);
    % self.colorbar_max=str2double(self.colorbar_max_string);
    
    % change the colorbar
    cb_min=self.model.colorbar_min;
    cb_max=self.model.colorbar_max;
    set(self.colorbar_axes_h,'YLim',[cb_min cb_max]);
    % cb_increment=(cb_max-cb_min)/256;
    % set(self.colorbar_h,'YData',[cb_min+0.5*cb_increment...
    %                              cb_max-0.5*cb_increment]);
    set(self.colorbar_h,'YData',[cb_min cb_max]);
    
    % % reset the z_slice index
    % self.z_slice_index=1;
    
    % prepare the axes to hold the z_slice
    is_a_file_open = self.model.is_a_file_open ;
    if is_a_file_open ,
        indexed_z_slice = self.model.indexed_z_slice ;
        [n_row,n_col]=size(indexed_z_slice);
        set(self.image_axes_h,'XLim',[0.5,n_col+0.5],...
            'YLim',[0.5,n_row+0.5]);
        
        % set the image to the current z_slice
        if ~isempty(self.image_h) ,
            delete(self.image_h)
        end
        self.image_h = ...
            image('Parent',self.image_axes_h,...
            'CData',indexed_z_slice,...
            'SelectionHighlight','off',...
            'ButtonDownFcn',@(~,~)(self.mouse_button_down_in_image()));
    else
        if ~isempty(self.image_h) ,
            delete(self.image_h)
        end
        self.image_h = [] ;
    end
    
    % % Need to update the roi border lines, labels to match the ones
    % % stored in the model!!
    % self.rois_changed_in_model();
    
    % update the z_slice counter stuff
    FPS_edit_string=sprintf('%g',self.model.fps) ;
    set(self.FPS_edit_h,'String',FPS_edit_string);
    
    if is_a_file_open ,
        set(self.z_index_edit_h, 'String', sprintf('%d',self.model.z_index));
        n_z_slice = self.model.z_slice_count ;
        set(self.of_z_slice_count_text_h, 'String', sprintf(' of %d',n_z_slice));
    else
        set(self.z_index_edit_h,'String', '');
        set(self.of_z_slice_count_text_h, 'String', '');
    end
    %of_z_slice_count_text_extent=get(self.of_z_slice_count_text_h,'Extent');
    %of_z_slice_count_text_position=get(self.of_z_slice_count_text_h,'Position');
    %pos_new=[of_z_slice_count_text_position(1:2) of_z_slice_count_text_extent(3:4)];
    %set(self.of_z_slice_count_text_h,'Position',pos_new);
    
    % Compute some things that effect what's enabled and what's not
    %at_least_one_roi_exists=~isempty(self.model.roi);
    
    % need to set image erase mode to none, since now there are no
    % more lines in front of the image
    %set(self.delete_all_rois_menu_h,'Enable',on_iff(at_least_one_roi_exists));
    %set(self.save_rois_to_file_menu_h,'Enable',on_iff(at_least_one_roi_exists));
    %set(self.hide_rois_menu_h,'Label',fif(self.hide_rois,'Show ROIs','Hide ROIs'));
    %set(self.hide_rois_menu_h,'Enable',on_iff(at_least_one_roi_exists));
    % disable the select ROI button
    %set(self.select_button_h,'Enable',on_iff(at_least_one_roi_exists));
    % enable the move all ROIs button
    %set(self.move_all_button_h,'Enable',on_iff(at_least_one_roi_exists));
    
    % enable buttons, menus that need enabling
    set(self.z_index_edit_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.FPS_edit_h,'enable',lapwing.on_iff(is_a_file_open));
    %set(self.elliptic_roi_button_h,'enable',lapwing.on_iff(a_video_is_open));
    %set(self.rect_roi_button_h,'enable',lapwing.on_iff(a_video_is_open));
    %set(self.polygonal_roi_button_h,'enable',lapwing.on_iff(a_video_is_open));
    set(self.zoom_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.pixel_data_type_min_max_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.min_max_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.five_95_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.abs_max_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.ninety_symmetric_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    %set(self.rois_menu_h,'enable',lapwing.on_iff(a_video_is_open));
    set(self.open_stack_menu_h,'enable',lapwing.on_iff(~is_a_file_open));
    set(self.close_stack_menu_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.to_start_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.play_backward_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.z_slice_backward_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.stop_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.z_slice_forward_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.play_forward_button_h,'enable',lapwing.on_iff(is_a_file_open));
    set(self.to_end_button_h,'enable',lapwing.on_iff(is_a_file_open));
    %set(self.mutation_menu_h,'enable',lapwing.on_iff(a_video_is_open));
    %set(self.motion_correct_menu_h,'enable',lapwing.on_iff(a_video_is_open));
    %set(self.load_overlay_menu_h,'enable',lapwing.on_iff(a_video_is_open));
    
    % % set the mode to elliptic_roi
    % self.set_mode('elliptic_roi');
    
end
