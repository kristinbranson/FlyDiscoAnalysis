function play(self, direction)

    % play the movie
    start_z_index = self.model.z_index ;
    z_slice_count = self.model.z_slice_count ;
    %n_rois=length(self.model.roi);
    % tempargh set(self.image_h,'EraseMode','none');
    fps = self.model.fps ;
    spf = 1/fps ;
    % if (direction>0)
    %   z_slice_sequence=start_z_index:z_slice_count;
    % else
    %   z_slice_sequence=start_z_index:-1:1;
    % end
    self.model.was_stop_button_hit = false ;
    z_index = start_z_index ;
    %for z_index=z_slice_sequence
    %tic;
    while (1<=z_index) && (z_index<=z_slice_count) ,
        %dt_this=toc;
        %fs=1/dt_this
        tic;
        self.model.z_index = z_index ;
        set(self.image_h, 'CData', self.model.indexed_z_slice) ;
        set(self.z_index_edit_h, 'String', num2str(z_index)) ;
        %self.sync_overlay();
        drawnow ;  % N.B.: this allows other callbacks to run!
        while (toc < spf)
        end
        if self.model.was_stop_button_hit ,
            break
        end
        z_index = z_index+direction ;
    end
    self.model.was_stop_button_hit = false ;

end
