%% set up path


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

%% get all experiments

data = SAGEListBowlExperiments('checkflags',false,'rootdir',rootdatadir,'daterange',{'20110201T000000'});


%%

successes = false(1,numel(data));
msgs = cell(1,numel(data));

%%

matlabpool(8)

%%

% disable waitbar
global DISABLESAGEWAITBAR;
DISABLESAGEWAITBAR = true; %#ok<NASGU>

%fid = fopen('AutomaticChecksSageResults20110719.txt','w');
parfor i = 104:numel(data),
  %fprintf('%d: %s...\n',i,data(i).experiment_name);
  expdir = data(i).file_system_path;
  [successes(i),msgs{i}] = FlyBowlAutomaticChecks_Sage(expdir,'dosave',false);
  %s = sprintf('%s\\n',msgs{i}{:});
  %s = s(1:end-2);
  %s = regexprep(s,',','.');
  %successletter = 'F';
  %if successes(i),
  %  successletter = 'P';
  %end
  %fprintf(fid,'%s,%s,%s\n',data(i).experiment_name,successletter,s);
  
end

%fclose(fid);
clear global DISABLESAGEWAITBAR;

%% print to file

fid = fopen('AutomaticChecksSageResults20110719.txt','w');
for i = 1:numel(data),
  if datenum(data(i).exp_datetime,'yyyymmddTHHMMSS') > 734695,
    continue;
  end
  msgs{i} = unique(msgs{i});
  msgs{i}(strcmp(msgs{i},'NULL')) = [];
  s = sprintf('%s\\n',msgs{i}{:});
  s = s(1:end-2);
  s = regexprep(s,',','.');
  successletter = 'F';
  if successes(i),
    successletter = 'P';
  end
  fprintf(fid,'%s,%s,%s\n',data(i).experiment_name,successletter,s);
end
fclose(fid);


%% compare to old stuff

oldfilename = 'AutomaticChecksSageResults20110616.txt';
fid = fopen(oldfilename,'r');
old_experiment_names = {};
old_successes = false(1,0);
while true,
  s0 = fgetl(fid);
  if ~ischar(s0),
    break;
  end
  s1 = regexp(s0,',','split');
  if numel(s1) ~= 2,
    error('Error parsing line %s',s0);
  end
  old_experiment_names{end+1} = s1{1}; %#ok<SAGROW>
  old_successes(end+1) = s1{2} == 'P'; %#ok<SAGROW>
end
fclose(fid);

nmatches = 0;
nmismatches = 0;
ismismatch = [];
for i = 1:numel(data),
  j = find(strcmp(data(i).experiment_name,old_experiment_names),1);
  if isempty(j), continue; end
  if old_successes(j) == successes(i),
    nmatches = nmatches + 1;
  else
    nmismatches = nmismatches + 1;
    
    s = sprintf('%s\\n',msgs{i}{:});
    s = s(1:end-2);
    s = regexprep(s,',','.');
    successletter = 'F';
    if successes(i),
      successletter = 'P';
    end
    
    fprintf('%d: %s,%s,%s\n',i,data(i).experiment_name,successletter,s);
  end
end