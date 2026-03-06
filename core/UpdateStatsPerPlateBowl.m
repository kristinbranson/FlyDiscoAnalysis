%% get bkgd stats per bowl from a lot of data

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

fns = {'bkgd_diagnostics_mean_bkgdcenter','bkgd_diagnostics_mean_bkgdcenter_llr',...
  'registrationdata_circleCenterX','registrationdata_circleCenterY','registrationdata_bowlMarkerTheta',...
  'temperature_diagnostics_mean','temperature_diagnostics_max'};

%% pull data

data = SAGEGetBowlData('dataset','score','checkflags',true,'removemissingdata',false,...
  'screen_type','primary','data_type',fns);

%% find and store median for each plate-bowl

[platebowls_bkgd,~,platebowlidx_bkgd] = unique({data.plate_bowl});
nbins = 50;
MINNCOUNTS = 20;
colors_bkgd = jet(numel(platebowls_bkgd))*.8;
for j = 1:numel(fns),
  fn = fns{j};
  y_bkgd = nan(1,numel(data));
  goodidx = ~cellfun(@isempty,{data.(fn)});
  y_bkgd(goodidx) = [data(goodidx).(fn)];
  edges = linspace(min(y_bkgd),max(y_bkgd),nbins+1);
  centers = (edges(1:end-1)+edges(2:end))/2;
  median_platebowl_bkgd.(fn) = nan(1,numel(platebowls_bkgd));
  figure(j + 10);
  clf;
  hold on;
  didplot = false(1,numel(platebowls_bkgd));
  centers1 = unique(y_bkgd);
  if numel(centers1) > nbins,
    centers1 = centers;
  end
  for i = 1:numel(platebowls_bkgd),
    y_curr = y_bkgd(platebowlidx_bkgd == i & goodidx);
    median_platebowl_bkgd.(fn)(i) = median(y_curr);
    counts = hist(y_curr,centers1);
    Z = sum(counts);
    if Z < MINNCOUNTS,
      disp(Z);
      continue;
    end
    counts = counts / Z;
    plot(centers1,counts,'.-','color',colors_bkgd(i,:));
    didplot(i) = true;
  end
  title(fn,'interpreter','none');
  legend(platebowls_bkgd(didplot),'interpreter','none');
end
fid = fopen(fullfile(settingsdir,analysis_protocol,'StatsPerPlateBowl.txt'),'w');
fprintf(fid,'platebowls');
fprintf(fid,',%s',platebowls_bkgd{:});
for i = 1:numel(fns),
  fn = fns{i};
  fprintf(fid,'\n%s',fn);
  fprintf(fid,',%f',median_platebowl_bkgd.(fn));
end
fprintf(fid,'\n');
fclose(fid);
