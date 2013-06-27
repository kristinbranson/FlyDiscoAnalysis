% subsample the data

expdir0 = '/groups/branson/bransonlab/tracking_data/CourtshipBowl/20120719_shelby';
expdir1 = '/groups/branson/bransonlab/tracking_data/CourtshipBowl/20120719_shelby_33fps';

td0 = load(fullfile(expdir0,'trx.mat'));

% subsample every 3rd frame
td1 = td0;
firstframetrack = min([td0.trx.firstframe]);
lastframetrack = max([td0.trx.endframe]);
tsample = [fliplr(firstframetrack:-3:1),firstframetrack+3:3:numel(td0.timestamps)];
ttracksample = firstframetrack:3:lastframetrack;
td1.timestamps = td0.timestamps(tsample);
nframestrack = numel(ttracksample);
td1.firstframetrack = find(tsample==ttracksample(1),1);
td1.lastframetrack = find(tsample==ttracksample(end),1);

idx = 1:3:numel(td0.trx(1).x);
fns = {'x','y','a','b','theta','sex'};
for i = 1:numel(fns),
  for j = 1:numel(td0.trx),
    td1.trx(j).(fns{i}) = td0.trx(j).(fns{i})(idx);
  end
end
dt = diff(td1.timestamps(td1.firstframetrack:td1.lastframetrack));
for j = 1:numel(td0.trx),
  td1.trx(j).dt = dt;
  td1.trx(j).firstframe = td1.firstframetrack;
  td1.trx(j).endframe = td1.lastframetrack;
  td1.trx(j).nframes = nframestrack;
  td1.trx(j).off = 1 - td1.trx(j).firstframe;
  td1.trx(j).timestamp = td0.timestamps(td0.firstframetrack:3:td0.lastframetrack);
end

save(fullfile(expdir1,'trx.mat'),'-struct','td1');

%% do the same thing for the wing tracking file

td0 = load(fullfile(expdir0,'wingtrx.mat'));

% subsample every 3rd frame
td1 = td0;
firstframetrack = min([td0.trx.firstframe]);
lastframetrack = max([td0.trx.endframe]);
tsample = [fliplr(firstframetrack:-3:1),firstframetrack+3:3:numel(td0.timestamps)];
ttracksample = firstframetrack:3:lastframetrack;
td1.timestamps = td0.timestamps(tsample);
nframestrack = numel(ttracksample);
td1.firstframetrack = find(tsample==ttracksample(1),1);
td1.lastframetrack = find(tsample==ttracksample(end),1);

idx = 1:3:numel(td0.trx(1).x);
fns = {'x','y','a','b','theta','sex','wing_anglel','wing_angler','xwingl','ywingl','xwingr','ywingr'};
for i = 1:numel(fns),
  for j = 1:numel(td0.trx),
    td1.trx(j).(fns{i}) = td0.trx(j).(fns{i})(idx);
  end
end
dt = diff(td1.timestamps(td1.firstframetrack:td1.lastframetrack));
for j = 1:numel(td0.trx),
  td1.trx(j).dt = dt;
  td1.trx(j).firstframe = td1.firstframetrack;
  td1.trx(j).endframe = td1.lastframetrack;
  td1.trx(j).nframes = nframestrack;
  td1.trx(j).off = 1 - td1.trx(j).firstframe;
  td1.trx(j).timestamp = td0.timestamps(firstframetrack:3:lastframetrack);
end

save(fullfile(expdir1,'wingtrx.mat'),'-struct','td1');

% and the per-frame features
if ~exist(fullfile(expdir1,'perframe'),'dir'),
  mkdir(fullfile(expdir1,'perframe'));
end
fns = {'nwingsdetected.mat','wing_areal.mat','wing_arear.mat','wing_trough_angle.mat'};

for i = 1:numel(fns),
  pd0 = load(fullfile(expdir0,'perframe',fns{i}));
  pd1 = pd0;
  for j = 1:numel(pd0.data),
    pd1.data{j} = pd0.data{j}(idx);
  end
  save(fullfile(expdir1,'perframe',fns{i}),'-struct','pd1');
end

%% register

td0 = load(fullfile(expdir1,'trx.mat'));
td1 = td0;
pxpermm = 6.9930;
fns = {'x','y','a','b'};
for i = 1:numel(fns),
  fnmm = [fns{i},'_mm'];
  for j = 1:numel(td0.trx),
    td1.trx(j).(fnmm) = td0.trx(j).(fns{i})/pxpermm;
  end
end
fps = 1/median(diff(td0.timestamps));
for j = 1:numel(td0.trx),
  td1.trx(j).pxpermm = pxpermm;
  td1.trx(j).fps = fps;
  td1.trx(j).theta_mm = td0.trx(j).theta;
end
save(fullfile(expdir1,'registered_trx.mat'),'-struct','td1');


%% register wing data

td0 = load(fullfile(expdir1,'wingtrx.mat'));
td1 = td0;
pxpermm = 6.9930;
fns = {'x','y','a','b'};
for i = 1:numel(fns),
  fnmm = [fns{i},'_mm'];
  for j = 1:numel(td0.trx),
    td1.trx(j).(fnmm) = td0.trx(j).(fns{i})/pxpermm;
  end
end
fps = 1/median(diff(td0.timestamps));
for j = 1:numel(td0.trx),
  td1.trx(j).pxpermm = pxpermm;
  td1.trx(j).fps = fps;
  td1.trx(j).theta_mm = td0.trx(j).theta;
end
save(fullfile(expdir1,'wingtracking_results.mat'),'-struct','td1');
