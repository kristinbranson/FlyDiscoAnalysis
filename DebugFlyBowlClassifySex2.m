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
  %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/expdirs';
  %rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110407';
  expdir = fullfile(rootdir,'GMR_12E07_AE_01_TrpA_Rig2Plate14BowlA_20110209T134320');
end

%% parameters

analysis_protocol = '20110804';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 

FlyBowlClassifySex2(expdir,params{:});

%% 

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

clear summary_diagnostics;
pfemale = cell(1,numel(expdirs));
medianarea = cell(1,numel(expdirs));
sex = cell(1,numel(expdirs));
areasmooth = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  try
    [trx,summary_diagnostics(i),areasmooth{i}] = FlyBowlClassifySex2(expdirs{i},params{:},'dosave',true); %#ok<SAGROW>
    sex{i} = {trx.sex};
    pfemale{i} = nan(1,numel(trx));
    medianarea{i} = nan(1,numel(trx));
    for fly = 1:numel(trx),
      pfemale{i}(fly) = nnz(strcmpi(trx(fly).sex,'F')) / trx(fly).nframes;
      medianarea{i}(fly) = nanmedian(areasmooth{i}{fly});
    end
  catch ME,
    getReport(ME)
  end
end

%%

clf;
hax = createsubplots(3,1,.05);
axes(hax(1)); %#ok<*MAXES>
h = nan(1,2);
h(1) = plot(1:numel(expdirs),[summary_diagnostics.classifier_mu_area_female],'r.-');
hold on;
h(2) = plot(1:numel(expdirs),[summary_diagnostics.classifier_mu_area_male],'b.-');
x = cell(1,numel(expdirs));
s = cell(1,numel(expdirs));
for i = 1:numel(expdirs),
  x{i} = zeros(1,numel(medianarea{i}))+i;
  s{i} = cellfun(@numel,areasmooth{i});
end
scatter(cell2mat(x),cell2mat(medianarea),cell2mat(s)*50/max(cell2mat(s)),cell2mat(pfemale),'.');
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