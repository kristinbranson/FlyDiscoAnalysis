function err = ComputeSetError(set_name,linestats,setstats,allstats,metadata,varargin)

statfns = {
  'velmag_ctr_flyany_frameany'
  'velmag_ctr_flyany_framewalk'
  'fractime_flyany_framestop'
  'fractime_flyany_framewinggrooming'
  'fractime_flyany_framewalk'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framebackup'
  'absdtheta_flyany_frameany'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framejump'
  'fractime_flyany_framerighting'
  'dnose2ell_flyany_frameany'
  'dcenter_flyany_frameany'
  'dcenter_flyany_framewalk'
  'dcenter_flyany_framestop'
  'nflies_close_flyany_frameany'
  'fractime_flyany_frametouch'
  'fractime_flyany_framechase'
  'fractime_flyany_frameattemptedcopulation'
  'wing_angle_diff_flyany_frameany'
  'fractime_flyany_framewingextension'
  'fractime_flyany_framewingflick'
  'dist2wall_flyany_frameany'
  'fractime_flyany_framenotanybehavior'
  };

[statfns,zplus] = myparse(varargin,'statfns',statfns,'zplus',[]);

nfns = numel(statfns);

if isempty(zplus),
  zplus = zeros(1,nfns);
end

err = nan(1,nfns);
seti = find(strcmp({setstats.metadata.set},set_name),1);
line_name = setstats.metadata(seti).line_name;
linei = find(strcmp(linestats.line_names,line_name));
maxnsets = linestats.nsets.velmag_ctr_flyany_frameany(linei);

expidx = strcmp({metadata.line_name},line_name) & ~strcmp({metadata.set},set_name);

if maxnsets <= 1,
  return;
end

  
for fni = 1:nfns,

  fn = statfns{fni};
  
  nsets = linestats.nsets.(fn)(linei);
  if nsets <= 1,
    continue;
  end
  setx = setstats.means.(fn)(seti);
  linex = (linestats.means.(fn)(linei)*nsets - setx)/(nsets-1);

  expstd = nanstd(allstats.(fn)(expidx),1);
  if isnan(expstd) || expstd <= 0,
    expstd = 1;
  end
  
  err(fni) = abs(linex-setx)/(expstd+zplus(fni));
  
end
