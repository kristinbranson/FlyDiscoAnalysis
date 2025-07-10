function point=nearest_visible_corner(self,current_point)

% get just the x and y
point=[current_point(1,1) ; current_point(1,2)];

% change to 'corner' coords
point=round(point-0.5)+0.5;

% get to the visible image dimensions
image_axes_h=findobj(self.figure_h,'Tag','image_axes_h');
xlim=get(image_axes_h,'XLim');
ylim=get(image_axes_h,'YLim');

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
