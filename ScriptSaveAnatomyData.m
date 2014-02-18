%% create anndata

% load in image choice data
load('/groups/branson/bransonlab/projects/olympiad/cross_assay/trunk/matlab/kristin/ImageChoices20130608.mat');
[line_names,~,lineidx] = unique({imdata.line});

% load in mask
maskfile = 'FullBrainMaskSymmetric.mat';
load(maskfile,'mask');
mask = permute(mask,[2,1,3]);

save AnnData20140204.mat imagechoice imdata line_names mask;

anndata = load('AnnData20140204.mat');
datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat';
load(datafile,'line_names');

%% data for params file

savedir = '/groups/branson/bransonlab/projects/olympiad/AverageAnatomyData20140204';
useallims = true;
maxqi = .1;
filsizexy = 5;
filsizez = 2;
cm = kjetsmooth(256);
save(fullfile(savedir,'params.mat'),'savedir','useallims','maxqi','filsizexy','filsizez','cm','anndata');

%% apply to lines

for i = 1:numel(line_names_curr),
  if ~ismember(line_names_curr{i},anndata.line_names), continue; end
  meanim = AverageLineAnatomy(line_names_curr{i},fullfile(anatomydatadir,'params.mat'));
end

%% apply to the rest of the lines

for i = 1:numel(line_names),
  filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
  if exist(filename,'file'),
    continue;
  end
  if ~ismember(line_names{i},anndata.line_names), 
    continue; 
  end
  meanim = AverageLineAnatomy(line_names{i},fullfile(anatomydatadir,'params.mat'));
end

%% compute mean image over all lines

nlinesread = 0;
for i = 1:numel(line_names),
  
  fprintf('Processing data for line %d / %d...\n',i,numel(line_names));
  
  filename = fullfile(savedir,sprintf('meanim_%s.mat',line_names{i}));
    
  if ~exist(filename,'file'),
    warning('File %s does not exist.',filename);
    continue;
  end
  datacurr = load(filename,'meanim');
  if nlinesread == 0,
    meanim = double(datacurr.meanim);
  else
    meanim = meanim + double(datacurr.meanim);
  end
  nlinesread = nlinesread+1;  
  
end

meanim = meanim / nlinesread;

mkdir(fullfile(savedir,'alllines_anatomy_20140209'));
save(fullfile(savedir,'alllines_anatomy_20140209','meanim.mat'),'meanim','line_names');