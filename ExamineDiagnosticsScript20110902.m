%% set up paths
[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    addpath E:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('E:\Code\cross_assay\matlab'));
    addpath(genpath('E:\Code\box\PostAnalysis'));
    addpath E:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',

    addpath C:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Code\cross_assay\trunk\matlab'));
    addpath(genpath('C:\Code\box\PostAnalysis\trunk\matlab'));
    addpath C:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';

  case 'bransonk-desktop',
    
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
  case 'robiea-ww1'
    
    addpath C:\Users\robiea\Documents\Code_versioned\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\cross_assay'));
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\PostAnalysis'));
    addpath C:\Users\robiea\Documents\Code_versioned\JCtrax\misc;
    addpath C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\settings';
    
  case 'robiea-ws'
    addpath '/groups/branson/home/robiea/Code_versioned/SAGE/MATLABInterface/Trunk';
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/matlab'));
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/postanalysis'));
    rmSvnPath
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/misc';
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/filehandling';
    addpath '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/olydat_browser';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/settings';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
end

analysis_protocol = '20110804';

%% pull some data


params = struct;
% data set: data -- we want behavior data
params.dataset = 'score';
% only pull data marked as primary screen data
params.screen_type = 'primary';
% don't use the automatic flag checking implemented in SAGEGetBowlData,
% we'll do that here if we want to
params.checkflags = false;
% don't remove experiments with missing data, we'll do that manually if we
% decide to. 
params.removemissingdata = false;
% no constraint on daterange, but you would implement it like this if you
% wanted one
params.daterange = {'20110824T000000','20120101T000000'};
% don't pull aborted experiments
params.flag_aborted = 0;
% don't pull experiments with flag_redo set to 1
params.flag_redo = 0;
% restrict to experiments with automated_pf = P
% params.automated_pf = 'P';
% restrict to experiments with manual_pf ~= F
% params.manual_pf = {'P','U'};

paramscell = struct2paramscell(params);

if false,
  data = SAGEGetBowlData(paramscell{:});
else
  load ExamineDiagnosticsData20110905.mat data;
end

x = datenum({data.exp_datetime},'yyyymmddTHHMMSS');
goodidx = strcmpi({data.automated_pf},'P') & ~strcmpi({data.manual_pf},'F');

for i = 1:numel(data),
  data(i).rig_plate_bowl = sprintf('%d.%s%d',data(i).rig,data(i).bowl,data(i).plate);
end

[rigplatebowls,~,rigplatebowlidx] = unique({data.rig_plate_bowl});
colors = jet(numel(rigplatebowls))*.8;

fprintf('failed experiments:\n');
fprintf('%s\n',data(~goodidx).experiment_name);

%% examine background diagnostics: brightness of image

hfig = 1;
figure(hfig);

% look at mean_bkgdcenter -- this is like the brightness of the image
clf;
subplot(2,1,1);
hold on;
y = nan(1,numel(data));
y(goodidx) = [data(goodidx).bkgd_diagnostics_mean_bkgdcenter];
for i = 1:numel(rigplatebowls),
  idx = rigplatebowlidx == i & goodidx;
  plot(x(idx),y(idx),'.','color',colors(i,:));
end
datetick('x',26);
title('Mean bkgd pixel intensity');
ylabel('Pixel intensity');

% look at difference from the median for that bowl
subplot(2,1,2);
hold on;
for i = 1:numel(rigplatebowls),
  idx = rigplatebowlidx == i & goodidx;
  y(idx) = y(idx) - mean(y(idx));
  plot(x(idx),y(idx),'.','color',colors(i,:));
end
datetick('x',26);
legend(rigplatebowls);
title('Normalized mean bkgd pixel intensity');
ylabel('Pixel intensity difference');

%% get bkgd stats per bowl from a lot of data

fns = {'bkgd_diagnostics_mean_bkgdcenter','bkgd_diagnostics_mean_bkgdcenter_llr'};
bkgddata = SAGEGetBowlData('dataset','score','checkflags',true,'removemissingdata',true,...
  'screen_type','primary','data_type',fns);
%% find and store median for each plate-bowl
[platebowls_bkgd,~,platebowlidx_bkgd] = unique({bkgddata.plate_bowl});
nbins = 100;
colors_bkgd = jet(numel(platebowls_bkgd))*.8;
for j = 1:numel(fns),
  fn = fns{j};
  y_bkgd = nan(1,numel(bkgddata));
  goodidx = ~cellfun(@isempty,{bkgddata.(fn)});
  y_bkgd(goodidx) = [bkgddata(goodidx).(fn)];
  edges = linspace(min(y_bkgd),max(y_bkgd),nbins+1);
  centers = (edges(1:end-1)+edges(2:end))/2;
  median_platebowl_bkgd.(fn) = nan(1,numel(platebowls_bkgd));
  figure(j + 10);
  clf;
  hold on;
  for i = 1:numel(platebowls_bkgd),
    y_curr = y_bkgd(platebowlidx_bkgd == i & goodidx);
    median_platebowl_bkgd.(fn)(i) = median(y_curr);
    counts = hist(y_curr,centers);
    counts = counts / sum(counts);
    plot(centers,counts,'.-','color',colors_bkgd(i,:));
  end
  title(fn,'interpreter','none');
  legend(platebowls_bkgd,'interpreter','none');
end
fid = fopen(fullfile(settingsdir,analysis_protocol,'BkgdStatsPerPlateBowl.txt'),'w');
fprintf(fid,'platebowls');
fprintf(fid,',%s',platebowls_bkgd{:});
for i = 1:numel(fns),
  fn = fns{i};
  fprintf(fid,'\n%s',fn);
  fprintf(fid,',%f',median_platebowl_bkgd.(fn));
end
fprintf(fid,'\n');
fclose(fid);


%% look at background ll

hfig = 2;
figure(hfig);

% look at bkgd ll -- how much the image looks like the experiment-wide
% model of background
clf;
subplot(2,1,1);
hold on;
y = nan(1,numel(data));
y(goodidx) = [data(goodidx).bkgd_diagnostics_mean_bkgdcenter_llr];
for i = 1:numel(rigplatebowls),
  idx = rigplatebowlidx == i & goodidx;
  plot(x(idx),y(idx),'.','color',colors(i,:));
end
datetick('x',26);
title('Bkgd ll');
ylabel('Bkgd ll');

% look at difference from the median for that bowl
subplot(2,1,2);
hold on;
for i = 1:numel(rigplatebowls),
  idx = rigplatebowlidx == i & goodidx;
  y(idx) = y(idx) - mean(y(idx));
  plot(x(idx),y(idx),'.','color',colors(i,:));
end
datetick('x',26);
legend(rigplatebowls);
title('Normalized bkgd ll');
ylabel('Normalized bkgd ll');

