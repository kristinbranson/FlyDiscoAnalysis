function update_pointer(self)

persistent above_image;

% init above_image
if isempty(above_image)
  mode=self.model.mode;
  if strcmp(mode,'select') || strcmp(mode,'move_all')
    above_image=false;
    set(self.figure_h,'Pointer','arrow');
  else  
    image_axes_h=self.image_axes_h;
    cp=get(image_axes_h,'CurrentPoint');
    cp=cp(1,1:2);
    xlim=get(image_axes_h,'XLim');
    ylim=get(image_axes_h,'YLim');
    above_image=~isempty(self.image_h) && ...
                xlim(1)<=cp(1) && ...
                cp(1)<=xlim(2) && ...
                ylim(1)<=cp(2) && ...
                cp(2)<=ylim(2);
    if above_image
      set(self.figure_h,'Pointer','crosshair');
    else
      set(self.figure_h,'Pointer','arrow');
    end
  end
end

% if we are above the image now and weren't before, or vice-versa, change
% the pointer appropriately
mode='select';
if ~( strcmp(mode,'select') || strcmp(mode,'move_all') )
  image_axes_h=self.image_axes_h;
  cp=get(image_axes_h,'CurrentPoint');
  cp=cp(1,1:2);
  xlim=get(image_axes_h,'XLim');
  ylim=get(image_axes_h,'YLim');
  above_image_now=~isempty(self.image_h) && ...
                  xlim(1)<=cp(1) && ...
                  cp(1)<=xlim(2) && ...
                  ylim(1)<=cp(2) && ...
                  cp(2)<=ylim(2);
  if above_image_now~=above_image
    if above_image_now
      set(self.figure_h,'Pointer','crosshair');
    else
      set(self.figure_h,'Pointer','arrow');
    end
    above_image=above_image_now;
  end
end

