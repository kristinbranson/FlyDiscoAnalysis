[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',
    
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-desktop',

    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

end

%% parameters

analysis_protocol = '20130722';
dateformat = 'yyyymmddTHHMMSS';
% params = {'settingsdir',settingsdir,...
%   'analysis_protocol',analysis_protocol,...
%   'daterange',{},...
%   'experiment_protocol',{'*8*','*7*'},...
%   'rearing_protocol','*8*'};

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'screen_type','primary',...
  'manual_pf',{'U','P'},...
  'automated_pf',{'U','P'},...
  'flag_redo','0',...
  'flag_aborted','0',...
  'checkflags',false};

datalocparamsfilestr = 'dataloc_params.txt';
startyear = 2011;
startmonth = 3;

%% load in metadata from file

load('CollectedPrimaryMetadata20130613.mat','metadata');
usemetadata = true;
idx = strcmp({metadata.line_name},'pBDPGAL4U');
metadata = metadata(idx);

%% loop through months


datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
timestamp = datestr(now,30);


year = startyear;
month = startmonth;

if usemetadata,
  dns = datenum({metadata.exp_datetime},dateformat);
end


while true,
  dvstart = [year,month,1,0,0,0];
  if month == 12,
    month = 1;
    year = year + 1;
  else
    month = month + 1;
  end
  dvend = [year,month,1,0,0,0];
  if datenum(dvend) > now,
    break;
  end

  daterange = {datestr(dvstart,dateformat),datestr(dvend,dateformat)};
  fprintf('Computing control data statistics for %s to %s...\n',daterange{:});
  outdir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,sprintf('%sto%s_computed%s',daterange{:},timestamp));
  
  if usemetadata,
    dnstart = datenum(dvstart);
    dnend = datenum(dvend);
    idx = dns>=dnstart & dns<=dnend;
    UpdateControlDataStats(rootdir,params{:},'daterange',daterange,'outdir',outdir,'experiments',metadata(idx));
  else
    UpdateControlDataStats(rootdir,params{:},'daterange',daterange,'outdir',outdir);
  end
end

%% make current link here

year = startyear;
month = startmonth;
pwdprev = pwd;
cd(dataloc_params.pBDPGAL4Ustatsdir);
while true,
  dvstart = [year,month,1,0,0,0];
  if month == 12,
    month = 1;
    year = year + 1;
  else
    month = month + 1;
  end
  dvend = [year,month,1,0,0,0];
  if datenum(dvend) > now,
    break;
  end
  daterange = {datestr(dvstart,dateformat),datestr(dvend,dateformat)};
  outdir = sprintf('%sto%s_computed%s',daterange{:},timestamp);
  if ~exist(outdir,'dir'),
    warning('%s does not exist',outdir);
    continue;
  end
  currentdir = sprintf('%sto%s_current',daterange{:});
  if exist(currentdir,'dir'),
    try
      unix(sprintf('rm %s',currentdir));
    catch ME,
      warning('Could not delete old current directory %s: %s',currentdir,getReport(ME));
    end
  end
  cmd = sprintf('ln -s %s %s',outdir,currentdir);
  unix(cmd);
end
cd(pwdprev);

%% also do everything into one directory

dvstart = [2011,02,1,0,0,0];
dvend = [2013,07,1,0,0,0];
daterange = {datestr(dvstart,dateformat),datestr(dvend,dateformat)};
fprintf('Computing control data statistics for %s to %s...\n',daterange{:});
outdir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,sprintf('%sto%s_computed%s',daterange{:},timestamp));
if usemetadata,
  UpdateControlDataStats(rootdir,params{:},'daterange',daterange,'outdir',outdir,'experiments',metadata);
else
  UpdateControlDataStats(rootdir,params{:},'daterange',daterange,'outdir',outdir);
end


%%

outdir = UpdateControlDataStats(rootdir,params{:});