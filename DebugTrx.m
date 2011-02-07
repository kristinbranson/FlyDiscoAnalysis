%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% data locations

if ispc,
  
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  analysis_protocol = '20110202_pc';
  dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,'dataloc_params.txt'));
  expdir = 'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734';
  expdir_read = fullfile(dataloc_params.rootreaddir,expdir);
  
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  protocol = 'CtraxTest20110202';
  analysis_protocol = '20110202';
  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol',protocol,'subreadfiles',{'ctrax_results.mat'});
  expdir = expdirs{1};
  expdir_read = expdir_reads{1};
end
% [readframe,nframes,fid,vidinfo] = get_readframe_fcn(fullfile(expdir_read,'movie.ufmf'));
% fclose(fid);

%% create Trx instance with one experiment
obj = Trx('settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
obj.AddExpDir(expdir);

%% try out all the per-frame properties

tmp = dir('compute_*.m');
fns = cellfun(@(s) s(9:end-2),{tmp.name},'UniformOutput',false);
fns_closestfly = {'magveldiff','veltoward','absthetadiff','absphidiff',...
  'anglefrom1to2','absanglefrom1to2'};
closestfly_types = {'center','nose2ell','ell2nose','anglesub'};
fns = setdiff(fns,fns_closestfly);
for i = 1:numel(fns_closestfly),
  for j = 1:numel(closestfly_types),
    fns{end+1} = sprintf('%s_%s',fns_closestfly{i},closestfly_types{j}); %#ok<SAGROW>
  end
end

for i = 1:numel(fns),
  fprintf('Computing %s\n',fns{i});
  plot(obj(1).(fns{i}));
  title(sprintf('%s ( %s )',fns{i},printunits(obj.units.(fns{i}))),'Interpreter','none');
  pause(2);
end

%% test to make sure things look okay

hfig = 1;
figure(hfig);
clf(hfig);
fly = 1;
ts = [100,3600,ceil(rand(1,2)*obj(fly).nframes)];
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'corfrac_min','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 2;
figure(hfig);
clf(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'angle2wall','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 3;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'dcenter','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 4;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'dnose2ell','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 5;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'dell2nose','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 6;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'anglesub','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 7;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'phi','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 8;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'magveldiff','hax',hax(i),'fly',fly,'t',ts(i));
end

hfig = 9;
figure(hfig);
fly = obj.nflies;
ts = ceil(rand(1,4)*obj(fly).nframes);
hax = createsubplots(1,numel(ts),.05,hfig);
for i = 1:numel(ts),
  IllustratePerFrameProps(obj,'veltoward','hax',hax(i),'fly',fly,'t',ts(i));
end