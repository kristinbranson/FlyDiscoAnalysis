function StoreTrajectories(obj,n,traj,dooverwrite)

if ~exist('dooverwrite','var'),
  dooverwrite = true;
end

traj_fns = {'x','y','theta','a','b','timestamps','area',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt','sex','area_mm',...
  'wing_anglel','wing_angler',...
  'xspine','yspine','xspine_mm','yspine_mm'};

for i = 1:numel(traj_fns),
  fn = traj_fns{i};
  if ~isfield(traj,fn), continue; end
  obj.SetPerFrameData(fn,{traj.(fn)},n);
  
  filename = obj.GetPerFrameFile(fn,n);
  data = {traj.(fn)}; %#ok<*NASGU>
  switch fn,
    case {'x','y','a','b'},
      units = parseunits('px');
    case {'x_mm','y_mm','a_mm','b_mm'},
      units = parseunits('mm');
    case {'theta','theta_mm','wing_anglel','wing_angler'},
      units = parseunits('rad');
    case 'dt',
      units = parseunits('s');
    case 'timestamps',
      units = parseunits('days');
    case 'sex',
      units = parseunits('unit');
  end
  try
    if dooverwrite && exist(filename,'file') && ~obj.DEBUG,
      try
        delete(filename);
      catch ME,
        warning('Could not delete file %s: %s',filename,getReport(ME));
      end
    end
    if (dooverwrite || ~exist(filename,'file')) && ~obj.DEBUG,
      save(filename,'data','units');
    end
  catch %#ok<CTCH>
    if ~exist(filename,'file'),
      error('Could not save to file %s',filename);
    else
      fprintf('Could not save to file %s, but file exists, so skipping\n',filename);
    end
  end
end
% conversion from pixels to mm
if isfield(traj,'pxpermm'),
  obj.pxpermm(n) = traj(1).pxpermm;
else
  obj.pxpermm(n) = nan;
end
% name of movie
if isfield(traj,'moviename'),
  obj.movienames{n} = traj(1).moviename;
end
flies = obj.exp2flies{n};
obj.firstframes(flies) = [traj.firstframe];
obj.endframes(flies) = [traj.endframe];
obj.nframes(flies) = [traj.nframes];
% fps
if isfield(traj,'fps'),
  obj.fps(n) = traj(1).fps;
else
  obj.fps(n) = nan;
end
