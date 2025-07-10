% check for video glitches

%% set up paths

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% parameters

datafile = 'ExperimentData_FlyBowl_20110618to20110625.mat';
nexps = 50;
nframes = 20;
suppradius = 25;
nax_r = 4;
nax_c = ceil(nframes/nax_r);

%% load data

load(datafile,'rawdata');
for i = 1:numel(rawdata),
  [~,expdir] = fileparts(rawdata(i).file_system_path);
  rawdata(i).file_system_path = fullfile(rootdatadir,expdir); %#ok<SAGROW>
end

%% sort by maxNBoxes

[~,order] = sort([rawdata.ufmf_diagnostics_summary_maxNBoxes],2,'descend');

%%

nboxes_all = cell(1,nexps);

%%
for ii = 1:nexps,
  
  disp(ii);
  i = order(ii);
  expdir = rawdata(i).file_system_path;
  headerinfo = ufmf_read_header(fullfile(expdir,'movie.ufmf'));
  nboxes = ufmf_read_nboxes(headerinfo,1:headerinfo.nframes);
  nboxes_all{ii} = nboxes;
  fclose(headerinfo.fid);
  
end

expdirs = {rawdata(order(1:nexps)).file_system_path};

save CheckForVideoGlitches.mat expdirs nboxes;

%% show results

hfig = 1;
figure(hfig);
clf;

for ii = 1:nexps,
  i = order(ii);
  expdir = rawdata(i).file_system_path;
  nboxes = nboxes_all{ii};
  
  moviename = fullfile(expdir,'movie.ufmf');
  while true,
    
    t = argmax(nboxes);
    uiwait(showufmf('UFMFName',moviename,'FirstFrame',t));
    nboxes(max(1,t-suppradius):min(numel(nboxes),t+suppradius)) = 0;
    res = questdlg('Show next frame or next movie or quit?','Next?','Frame','Movie','Quit','Frame');
    if strcmpi(res,'Movie'),
      break;
    end
    
  end
  
  if strcmpi(res,'Quit'),
    break;
  end
  
end