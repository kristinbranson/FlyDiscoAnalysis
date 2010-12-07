% heatmap = obj.CenterPositionHeatmap(varargin)

function heatmap = CenterPositionHeatmap(obj,varargin)

%% parse inputs

% parse optional inputs
[lim_x,lim_y,...
  conditions,...
  minspeed,...
  leftovers] = ...
  myparse_nocheck(varargin,...
  'lim_x',nan(1,2),'lim_y',nan(1,2),...
  'conditions',[],...
  'minspeed',0);

%% set lim_x, lim_y if not input

% left-most point in the arena
if isnan(lim_x(1)),
  lim_x(1) = obj.arena_center_mm(1)-obj.arena_radius_mm;
end
% right-most point in the arena
if isnan(lim_x(2)),
  lim_x(2) = obj.arena_center_mm(1)+obj.arena_radius_mm;
end
% top-most point in the arena
if isnan(lim_y(1)),
  lim_y(1) = obj.arena_center_mm(2)-obj.arena_radius_mm;
end
% bottom-most point in the arena
if isnan(lim_y(2)),
  lim_y(2) = obj.arena_center_mm(2)+obj.arena_radius_mm;
end

%% add minspeed to conditions
if minspeed > 0,
  if isempty(conditions),
    conditions = @(trk) [false,trk.velmag > minspeed];
  else
    conditions = @(trk) conditions(trk) & [false,trk.velmag > minspeed];
  end
end

%% use 2-D histogram to do everything else

heatmap = obj.HistogramTwoMeasurements('x_mm','y_mm',...
  'lim_x',lim_x,'lim_y',lim_y,'conditions',conditions,leftovers{:});

