function move_roi(self,action,roi_index_start)

persistent figure_h;
persistent image_axes_h;
persistent anchor;
persistent roi_border_h;
persistent roi_label_h;
persistent border_x_anchor;
persistent border_y_anchor;
persistent label_pos_anchor;
persistent roi_index;

switch(action)
  case 'start'
    figure_h=self.figure_h;
    image_axes_h=self.image_axes_h;
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    roi_border_hs=self.border_roi_h;
    roi_index=roi_index_start;
    roi_border_h=roi_border_hs(roi_index);
    border_x_anchor=get(roi_border_h,'XData');
    border_y_anchor=get(roi_border_h,'YData');
    roi_label_hs=self.label_roi_h;
    roi_label_h=roi_label_hs(roi_index);
    label_pos_anchor=get(roi_label_h,'Position');
    % set the callbacks for the drag
    set(figure_h,'WindowButtonMotionFcn',...
               @(~,~)(self.move_roi('move')));
    set(figure_h,'WindowButtonUpFcn',...
               @(~,~)(self.move_roi('stop')));
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    x=border_x_anchor+(point(1)-anchor(1));
    y=border_y_anchor+(point(2)-anchor(2));
    set(roi_border_h,'XData',x);
    set(roi_border_h,'YData',y);
    set(roi_label_h,'Position',label_pos_anchor+[point-anchor 0]);
  case 'stop'
    % change the move and buttonup calbacks
    set(figure_h,'WindowButtonMotionFcn', ...
                 @(~,~)(self.update_pointer()) );
    set(figure_h,'WindowButtonUpFcn',[]);
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    x=border_x_anchor+(point(1)-anchor(1));
    y=border_y_anchor+(point(2)-anchor(2));
    set(roi_border_h,'XData',x);
    set(roi_border_h,'YData',y);
    set(roi_label_h,'Position',label_pos_anchor+[point-anchor 0]);
    % store the new border in the model
    self.model.roi(roi_index).border=[x;y];
    % clear the persistents
    figure_h=[];
    image_axes_h=[];
    anchor=[];
    roi_border_h=[];
    roi_label_h=[];
    border_x_anchor=[];
    border_y_anchor=[];
    label_pos_anchor=[];
    roi_index=[];
end  % switch






