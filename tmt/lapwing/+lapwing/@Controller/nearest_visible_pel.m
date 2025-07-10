function point=nearest_visible_pel(self,current_point)

% get just the x and y
point=[current_point(1,1) ; current_point(1,2)];

% find nearest pel
point=round(point);

% get to the visible image dimensions
image_axes_h=findobj(self.figure_h,'Tag','image_axes_h');
xlim=get(image_axes_h,'XLim');
ylim=get(image_axes_h,'YLim');

% change xlim, ylim to pel center coords
xlim=[ceil(xlim(1)) floor(xlim(2))];
ylim=[ceil(ylim(1)) floor(ylim(2))];

% constrain the returned point to the visible image dims
if point(1)<xlim(1)    
  point(1)=xlim(1); 
end
if point(1)>xlim(2)
  point(1)=xlim(2); 
end
if point(2)<ylim(1)
  point(2)=ylim(1);
end
if point(2)>ylim(2)
  point(2)=ylim(2); 
end
