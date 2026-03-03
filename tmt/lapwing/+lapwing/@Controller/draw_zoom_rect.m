function draw_zoom_rect(self,action)

persistent figure_h;
persistent image_axes_h;
persistent anchor;
persistent rect_h;

% init the persistents that point to the image figure and the image 
% axes, if this is the first time through this function
if isempty(figure_h)
  figure_h=self.figure_h;
end  
if isempty(image_axes_h)
  image_axes_h=self.image_axes_h;
end  

switch(action)
  case 'start'
    point=self.nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    anchor=point;
    % create a new rectangle
    rect_h=...
      line('Parent',image_axes_h,...
           'Color',[0 0 1],...
           'Tag','border_h',...
           'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)],...
           'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)],...
           'ZData',[1 1 1 1 1],...
           'ButtonDownFcn','lapwing.callback');
    % set the callbacks for the drag
    set(figure_h,'WindowButtonMotionFcn',...
                 @(src,event)(self.draw_zoom_rect('move')));
    set(figure_h,'WindowButtonUpFcn',...
                 @(src,event)(self.draw_zoom_rect('stop')));
  case 'move'
    point=self.nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
  case 'stop'
    % change the move and buttonup calbacks
    set(figure_h,'WindowButtonMotionFcn', ...
                 @(src,event)(self.update_pointer()));
    set(figure_h,'WindowButtonUpFcn',[]);
    % now do the stuff we'd do for a move also
    point=self.nearest_visible_corner(get(image_axes_h,'CurrentPoint'));
    set(rect_h,...
        'XData',[anchor(1) anchor(1) point(1) point(1)  anchor(1)]);
    set(rect_h,...
        'YData',[anchor(2) point(2)  point(2) anchor(2) anchor(2)]);
    % clear the persistents
    figure_h=[];
    image_axes_h=[];
    % do the zoom
    delete(rect_h)
    self.zoom_in(point,anchor);
end  % switch






