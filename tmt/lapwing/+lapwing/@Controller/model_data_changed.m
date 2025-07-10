function model_data_changed(self)

% The data in the model (but not the dimensions of the data, nor the time
% information) has changed.  Need to update self appropriately.

% change the image data in the image HG object
set(self.image_h,'cdata',self.indexed_z_slice);

% update the figure
drawnow('update');
drawnow('expose');

end
