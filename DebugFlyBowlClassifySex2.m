%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath E:\Code\hmm;
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  expdir = 'E:\Data\FlyBowl\bowl_data\GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  expdir = fullfile(rootdir,'GMR_12E07_AE_01_TrpA_Rig2Plate14BowlA_20110209T134320');
end

%% parameters

analysis_protocol = '20110222';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 

FlyBowlClassifySex2(expdir,params{:});

%% 

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

clear summary_diagnostics;
pfemale = cell(1,numel(expdirs));
meanarea = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  try
    [trx,summary_diagnostics(i),areasmooth] = FlyBowlClassifySex2(expdirs{i},params{:},'dosave',true); %#ok<SAGROW>
    pfemale{i} = nan(1,numel(trx));
    meanarea{i} = nan(1,numel(trx));
    for fly = 1:numel(trx),
      pfemale{i}(fly) = nnz(strcmpi(trx(fly).sex,'F')) / trx(fly).nframes;
      meanarea{i}(fly) = nanmean(areasmooth{fly});
    end
  catch ME,
    getReport(ME)
  end
end

%%

clf;
hax = createsubplots(3,1,.05);
axes(hax(1));
h = nan(1,2);
h(1) = plot(1:numel(expdirs),[summary_diagnostics.classifier_mu_area_female],'r.-');
hold on;
h(2) = plot(1:numel(expdirs),[summary_diagnostics.classifier_mu_area_male],'b.-');
x = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  x{i} = 1:numel(meanarea{i});
end
scatter(cell2mat(x),cell2mat(meanarea),[],cell2mat(pfemale),'filled');
set(gca,'XTick',1:numel(expdirs),'XTickLabel',{});
legend(h,{'classifier_mu_area_female','classifier_mu_area_male'},...
  'interpreter','none','Location','best');
ylabel('Area (mm^2)');

axes(hax(2));
plot(1:numel(expdirs),[summary_diagnostics.mean_nfemales],'r.-');
hold on;
plot(1:numel(expdirs),[summary_diagnostics.mean_nmales],'b.-');
set(gca,'XTick',1:numel(expdirs),'XTickLabel',{});
ylabel('N. flies');
legend({'mean_nfemales','mean_nmales'},...
  'interpreter','none','Location','best');

axes(hax(3));
plot(1:numel(expdirs),[summary_diagnostics.mean_nswaps],'k.-');
ylabel('Mean n. swaps');

xticklabels = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  expinfo = parseExpDir(expdirs{i});
  if expinfo.line(1) == 'G',
    m = regexp(expinfo.line,'GMR_([^_]+)_.*','tokens','once');
    lineabbr = m{1};
  else
    lineabbr = 'pB';
  end
  xticklabels{i} = sprintf('%s %s%s',lineabbr,expinfo.rig,expinfo.bowl);
end
set(gca,'XTick',1:numel(expdirs),'XTickLabel',xticklabels);
xlabel('Experiment');