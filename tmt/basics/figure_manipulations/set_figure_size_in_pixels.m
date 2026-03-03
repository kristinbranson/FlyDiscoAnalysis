function set_figure_size_in_pixels(varargin)

if length(varargin)==1 ,
  f = gcf() ;
  sz_fig_new = varargin{1} ;
elseif length(varargin)==2 ,
  f = varargin{1} ;
  sz_fig_new = varargin{2} ;  
else
  error('Too many arguments') ;
end

% sz_fig_new in inches

% save screen, fig units
original_screen_units=get(0,'units');
original_fig_units=get(f,'units');

% set screen, figure units
set(0,'units','pixels');
set(gcf,'units','pixels');

% calculate new figure offset
pos_screen=get(0,'screensize');  sz_screen=pos_screen(3:4);
pos_fig=get(f,'position');
top_edge_fig=pos_fig(2)+pos_fig(4);
fig_offset_new=[(sz_screen(1)-sz_fig_new(1))/2 top_edge_fig-sz_fig_new(2)];

% set fig to new position
fig_pos_new=[fig_offset_new sz_fig_new];
set(f,'position',fig_pos_new);

% replace units
set(0,'units',original_screen_units);
set(f,'units',original_fig_units);
