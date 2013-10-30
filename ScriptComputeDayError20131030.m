alldates = unique(cellfun(@(x) x(1:8),{setstats.metadata.exp_datetime},'UniformOutput',false));

hfig = 1;
figure(hfig);
clf;
hold on;
hax = gca;
datetick('x');

statfnscheck = {
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

[~,idx] = ismember(statfnscheck,statfns);
zplus = controlstd(idx);

clear errsummary err;
for i = randperm(numel(alldates)),
  date = alldates{i};
  fprintf('Date %s (%d/%d)\n',date,i,numel(alldates));
  dn = datenum(date,'yyyymmdd');
  [errsummary(i),err{i}] = ComputeDayError(date,linestats,setstats,allstats,metadata,'statfns',statfnscheck,'zplus',zplus);
  plot(hax,dn,errsummary(i).maxerr,'k.');
  axisalmosttight(.01,hax);
  drawnow;
end
