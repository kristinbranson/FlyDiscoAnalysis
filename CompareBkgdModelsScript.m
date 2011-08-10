%% compare background models

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture/intensity';
end

%% experiments

expdirs = dir(fullfile(rootdatadir,'*20*'));
basenames = {expdirs.name};
expdirs = cellfun(@(s) fullfile(rootdatadir,s),basenames,'UniformOutput',false);

experiments = [];
for i = 1:numel(expdirs),
  experiments = structappend(experiments,parseExpDir(expdirs{i}));
end

[~,order] = sort({experiments.date});
expdirs = expdirs(order);
basenames = basenames(order);
experiments = experiments(order);

%% load quickstats

quickstats = [];
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  quickstats = structappend(quickstats,loadQuickStats(expdir));  
end

%% plot scan lines

figure(1);
clf;
hax = createsubplots(2,2,.05);

%colors = jet(numel(quickstats))*.75;
colors = [repmat([.7,0,0],[4,1]);zeros(4,3)];

for i = 1:4,
  fn = sprintf('BkgdScanLine_intensities_%d',i);
  axes(hax(i)); %#ok<LAXES>
  h = plot(cat(1,quickstats.(fn))','-');
  for j = 1:numel(quickstats),
    set(h(j),'color',colors(j,:));
  end

  axisalmosttight;
  if i == 1,
    for j = 1:4,
      legends{j} = sprintf('no cardboard %s%s',experiments(j).rig,experiments(j).bowl);
    end
    for j = 1:4,
      legends{4+j} = sprintf('with cardboard %s%s',experiments(4+j).rig,experiments(4+j).bowl);
    end
    legend(h,legends,'interpreter','none');
  end
  title(sprintf('Scan line %d',i));
end
