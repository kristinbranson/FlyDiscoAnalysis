function StoreTrajectories(obj,n,traj)

traj_fns = {'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','dt'};

for i = 1:numel(traj_fns),
  fn = traj_fns{i};
  if ~isfield(traj,fn), continue; end
  obj.SetPerFrameData(fn,{traj.(fn)},n);
  
  filename = obj.GetPerFrameFile(fn,n);
  data = traj.(fn); %#ok<*NASGU>
  switch fn,
    case {'x','y','a','b'},
      units = parseunits('px');
    case {'x_mm','y_mm','a_mm','b_mm'},
      units = parseunits('mm');
    case {'theta','theta_mm'},
      units = parseunits('rad');
    case 'dt',
      units = parseunits('s');
    case 'timestamps',
      units = parseunits('days');
  end
  save(filename,'data','units');
end
% conversion from pixels to mm
if isfield(traj,'pxpermm'),
  obj.pxpermm(n) = traj(1).pxpermm;
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
end
