expdir = '/groups/branson/bransonlab/taylora/flydisco/example-experiments/nice-early-pfly-experiment-with-tracking-new/pflySingleMale_pfly_SlowRamp_LexAChrimx71G01LexA_UASChrimsonVenusX0070_20240320T095847';
%%
td = load(fullfile(expdir,'registered_trx.mat'));
%%
fns_check0 = {'timestamps','x','y','theta','a','b','xwingl','ywingl','xwingr','ywingr','x_mm','y_mm','a_mm','b_mm','theta_mm'};
fns_check1 = {'dt'};
for i = 1:numel(td.trx),
  t = td.trx(i).nframes;
  for j = 1:numel(fns_check0),
    if numel(td.trx(i).(fns_check0{j})) ~= t,
      fprintf('%d: %s: %d != %d\n',i,fns_check0{j},numel(td.trx(i).fns_check0{j}),t)
    end
  end
  for j = 1:numel(fns_check1),
    if numel(td.trx(i).(fns_check1{j})) ~= t-1,
      fprintf('%d: %s: %d != %d\n',i,fns_check1{j},numel(td.trx(i).fns_check1{j}),t-1)
    end
  end
end
%%
pfns = dir(fullfile(expdir,'perframe','*.mat'));
for i = 1:numel(pfns),
  dcurr = load(fullfile(pfns(i).folder,pfns(i).name));
  for j = 1:numel(td.trx),
    t = td.trx(j).nframes;
    if size(dcurr.data{j},2) < t-1 || size(dcurr.data{j},2) > t,
      fprintf('%d: %s: %d != (%d,%d)\n',j,pfns(i).name,size(dcurr.data{j},2),t-1,t);
    end
  end
end