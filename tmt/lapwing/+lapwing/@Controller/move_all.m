function move_all(self,action)

persistent figure_h;
persistent image_axes_h;
persistent n_rois;
persistent roi_borders;
persistent roi_border_line_hs;
persistent roi_label_hs;
persistent anchor;

switch(action)
  case 'start'
    % set up the persistents
    figure_h=self.figure_h;
    image_axes_h=self.image_axes_h;
    roi_border_line_hs=self.border_roi_h;
    n_rois=length(roi_border_line_hs);
    roi_borders=cell(n_rois,1);
    for i=1:n_rois
      x=get(roi_border_line_hs(i),'XData');
      y=get(roi_border_line_hs(i),'YData');
      this_roi_border=zeros(2,length(x));
      this_roi_border(1,:)=x;
      this_roi_border(2,:)=y;
      roi_borders{i}=this_roi_border;
    end
    % hide the ROI labels
    roi_label_hs=self.label_roi_h;
    set(roi_label_hs,'Visible','off');
    % now do real stuff    
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2); 
    anchor=point;
    % set the callbacks for the drag
    set(figure_h,'WindowButtonMotionFcn',...
                 @(~,~)(self.move_all('move')) );
    set(figure_h,'WindowButtonUpFcn',...
                 @(~,~)(self.move_all('stop')) );
  case 'move'
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    vect=point-anchor;
    % move all the borders & labels
    for i=1:n_rois
      this_roi_border=roi_borders{i};
      set(roi_border_line_hs(i),'XData',vect(1)+this_roi_border(1,:));
      set(roi_border_line_hs(i),'YData',vect(2)+this_roi_border(2,:));
    end
  case 'stop'
    % change the move and buttonup calbacks
    set(figure_h,'WindowButtonMotionFcn',@(~,~)(self.update_pointer()));
    set(figure_h,'WindowButtonUpFcn',[]);
    % now do the stuff we'd do for a move also
    cp=get(image_axes_h,'CurrentPoint');
    point=cp(1,1:2);
    vect=point-anchor;
    % move all the borders & labels
    for i=1:n_rois
      this_roi_border=roi_borders{i};
      n_vertex=size(this_roi_border,2);
      r=repmat(vect',[1 n_vertex])+this_roi_border;
      set(roi_border_line_hs(i),'XData',r(1,:));
      set(roi_border_line_hs(i),'YData',r(2,:));
      % store the new border in the model
      self.model.roi(i).border=r;
    end
    % move the labels, make visible again
    for i=1:n_rois
      this_position=get(roi_label_hs(i),'Position');
      set(roi_label_hs(i),'Position',[vect 0]+this_position);
      set(roi_label_hs(i),'Visible','on');
    end        
    % clear the persistents
    figure_h=[];
    image_axes_h=[];
    n_rois=[];
    roi_borders=[];
    roi_border_line_hs=[];
    roi_label_hs=[];
    anchor=[];
    % since the borders are stored in the lines, we're done
end  % switch

