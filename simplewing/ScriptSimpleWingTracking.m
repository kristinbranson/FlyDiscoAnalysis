%% set up paths

rootdatadir = '/groups/branson/bransonlab/projects/JAABA/data/FlyBowl';
%rootdatadir = 'C:\Data\JAABA\FlyBowl';
addpath ../misc;
addpath ../filehandling;

%% parameters

% file names
experiment_name = 'GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818';
expdir = fullfile(rootdatadir,experiment_name);
paramsfile = 'WingTrackingParameters.txt';

DEBUG = 1;

%%

[trx,perframedata] = TrackWings(expdir,...
  'paramsfile',paramsfile,...
  'debug',DEBUG,...
  'outtrxfilestr','test.mat');

%% get a list of all experiments with Ctrax results

rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
requiredfiles = {'movie.ufmf'
  'movie.ufmf.ann'
  'registered_trx.mat'};

dirs = dir(rootdatadir);
fid = fopen('expdirs_trackwings20130603.txt','w');
expdirs = {};
for i = i:numel(dirs),
  if ~dirs(i).isdir,
    continue;
  end
  name = dirs(i).name;
  if ismember(name,{'.','..'}),
    continue;
  end
  metadatafile = fullfile(rootdatadir,name,'Metadata.xml');
  if ~exist(metadatafile,'file'),
    fprintf('No metadata for %s\n',name);
    continue;
  end
  try
    metadata = ReadMetadataFile(metadatafile);
  catch ME,
    fprintf('Could not read metadata from %s: %s\n',metadatafile,getReport(ME));
    continue;
  end
  if ~isfield(metadata,'screen_type') || ~strcmp(metadata.screen_type,'primary'),
    continue;
  end
  
  filesexist = true;
  for j = 1:numel(requiredfiles),
    if ~exist(fullfile(rootdatadir,name,requiredfiles{j}),'file'),
      missingfile = requiredfiles{j};
      filesexist = false;
      break;
    end
  end
  if ~filesexist,
    fprintf('Not all files exist for %s (missing %s)\n',name,missingfile);
    continue;
  end
  fprintf(fid,'%s\n',fullfile(rootdatadir,name));
  expdirs{end+1} = fullfile(rootdatadir,name);
  fprintf('Adding %d: %s\n',numel(expdirs),fullfile(rootdatadir,name));
end
fclose(fid);
