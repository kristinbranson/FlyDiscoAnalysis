% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootoutdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments';
analysis_protocol = '20120828_non_olympiad_robiea_age_20120519';
mintimestamp = '20120501T000000';
npar = 8;
datatype = 'wildtype2';


%% set data locations

if strcmpi(datatype,'starvation'),
  
  [exp_params,expdirs] = StarvationSAGEParams();
  outdatadir = fullfile(rootoutdatadir,'starvation');
  outexpdirfilename = 'expdirs_starvation20120909.txt';
  
elseif strcmpi(datatype,'age'),
  
  [exp_params,expdirs] = AgeSAGEParams();
  outdatadir = fullfile(rootoutdatadir,'age');
  outexpdirfilename = 'expdirs_age20120909.txt';
  
elseif strcmpi(datatype,'wildtype2'),
  
  [exp_params,expdirs] = WildtypeSAGEParams();
  outdatadir = fullfile(rootoutdatadir,'wildtype');
  outexpdirfilename = 'expdirs_wildtype20120914.txt';
  
end

metadata = [];
protocols = fieldnames(exp_params);
nprotocols = numel(protocols);
for i = 1:nprotocols,
  metadatacurr = SAGEListBowlExperiments(exp_params.(protocols{i}){:},'checkflags',false,'removemissingdata',false);
  [metadatacurr.condition] = deal(protocols{i});
  if i == 1,
    metadata = metadatacurr;
  else
    metadata = [metadata,metadatacurr]; %#ok<AGROW>
  end
end

badidx = [metadata.manual_pf] == 'F';
metadata(badidx) = [];

for i = 1:nprotocols,
  protocol = protocols{i};  
  idx = strcmp({metadata.condition},protocol);
  tmp1 = setdiff({metadata(idx).file_system_path},expdirs.(protocol));
  if ~isempty(tmp1),
    fprintf('%s: The following directories were pulled from SAGE but not listed manually\n',protocol);
    fprintf('%s\n',tmp1{:});
  end
  tmp2 = setdiff(expdirs.(protocol),{metadata(idx).file_system_path});
  if ~isempty(tmp2),
    fprintf('%s: The following directories were listed manually but not pulled from SAGE\n',protocol);
    fprintf('%s\n',tmp2{:});
  end
end

% copy over data
for i = 1:nprotocols,
  protocol = protocols{i};  
  idx = find(strcmp({metadata.condition},protocol));
  for j = idx(:)',
    path0 = metadata(j).file_system_path;
    [~,name] = fileparts(path0);
    path1 = fullfile(outdatadir,protocol,name);
    if exist(path1,'dir'),
      fprintf('%s already exists\n',path1);
      continue;
    end
    fprintf('copying %s to %s\n',path0,path1);
    SymbolicCopyExperimentDirectory(path0,fullfile(outdatadir,protocol));
  end
end

% update paths
for i = 1:nprotocols,
  protocol = protocols{i};
  idx = find(strcmp({metadata.condition},protocol));
  for j = idx(:)',
    path0 = metadata(j).file_system_path;
    [~,name] = fileparts(path0);
    path1 = fullfile(outdatadir,protocol,name);
    if ~exist(path1,'dir'),
      error('%s does not exist',path1);
    end
    metadata(j).file_system_path = path1; %#ok<SAGROW>
  end
end

% fid = fopen(fullfile(outdatadir,outexpdirfilename),'w');
% fprintf(fid,'%s\n',metadata.file_system_path);
% fclose(fid);


%%
% 
% metadata = SAGEListBowlExperiments('screen_type','non_olympiad_robiea*','checkflags',false,'removemissingdata',false);
% 
% for i = 1:numel(metadata),
%   fn = fullfile(metadata(i).file_system_path,'analysis_protocol.txt');
%   if ~exist(fn,'file'),
%     continue;
%   end
%   fid = fopen(fn,'r');
%   s = fgetl(fid);
%   s = strtrim(s);
%   fclose(fid);
%   metadata(i).analysis_protocol = s;
% end
% 
% % check analysis protocols
% [allprotocols,~,protocolidx] = unique({metadata.analysis_protocol});
% 
% % 282 experiments pulled. 
% % 44 had manual_pf = F
% % 1 had automated_pf = F for flies loaded time too long
% % analysis_protocol = NULL for all experiments with manual_pf = F
% % data collected in may to june has analysis_protocol
% % ctrax_version:0.3.1,ctrax_settings:20111221,analysis_settings:20120210,analysis_version:a407adaee6097e3ecb42a6645769b9c5fa234efd
% % data collected in august to september has analysis_protocol
% % ctrax_version:0.3.1,ctrax_settings:20111221,analysis_settings:20120706,analysis_version:c64e387fc87036c2cb66ab409a5865ff54233a7c
% 
% for i = 1:numel(metadata),
%   if metadata(i).automated_pf == 'F',
%     continue;
%   end
%   [~,name] = fileparts(metadata(i).file_system_path);
%   switch metadata(i).screen_type,
%     case 'non_olympiad_robiea_starvation',
%       cmd = sprintf('ls -d %s',fullfile('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation/starve*h',name));
%  
%     case 'non_olympiad_robiea_age',
%       cmd = sprintf('ls -d %s',fullfile('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age/age*days',name));
% 
%     case {'non_olympiad_robiea_wildtype','non_olympiad_robiea_wildtype2'},
%       cmd = sprintf('ls -d %s',fullfile('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/wildtype/results',name));
%     otherwise
%       error('unknown screen_type %s',metadata(i).screen_type);
%       
%   end
%   [tmp1,tmp] = unix(cmd);
%   tmp = strtrim(tmp);
%   if ~tmp1 && ~isempty(tmp),
%     metadata(i).file_system_path2 = tmp;
%   else
%     metadata(i).file_system_path2 = metadata(i).file_system_path;
%   end
%       
% end
% 
% badidx = [metadata.automated_pf] == 'F';
% metadata(badidx) = [];

%experiment_names = cellfun(@(x) x(9:end),{metadata([metadata.automated_pf]=='P').experiment_name},'UniformOutput',false);

%% initialize

issuccess = nan(1,numel(metadata));
msgs = cell(1,numel(metadata));
nproblems = 0;
ii = 0;
datenums = datenum({metadata.exp_datetime},'yyyymmddTHHMMSS');

hfig = 1;
figure(hfig);
clf;

hdata = plot(datenums,issuccess,'k.-');
set(gca,'XLim',[min(datenums)-1,max(datenums)+1],'YLim',[-.25,1.25]);
datetick('x','mmmyy','keeplimits');
hti = title(sprintf('%d problems / %d checked',nproblems,ii));

order = randperm(numel(metadata));

%% check

for ii = 1:npar:numel(metadata),
  
  %i = order(ii);
  
  idxcurr = order(ii:min(ii+npar-1,numel(metadata)));
  currissuccess = nan(1,numel(idxcurr));
  currmsgs = cell(1,numel(idxcurr));
  curr_expdirs = {metadata(idxcurr).file_system_path};
  
  parfor j = 1:numel(idxcurr),
    fprintf('***** Checking experiment %s...\n',curr_expdirs{j});
    [currissuccess(j),currmsgs{j}] = CheckExperiment20120726(curr_expdirs{j},'analysis_protocol',analysis_protocol);
  end
  
  issuccess(idxcurr) = currissuccess;
  msgs(idxcurr) = currmsgs;
  nproblems = nnz(issuccess==0);
  
  set(hdata,'YData',issuccess);
  set(hti,'String',sprintf('%d problems / %d checked',nproblems,ii+numel(idxcurr)-1));
  
  drawnow;
  
  if mod(ii,100) < npar,
    save(sprintf('CheckExperiments_robiea_%s_20120909.mat',datatype),'issuccess','msgs','order','ii','metadata');
  end
  
end

%% check for scores files

startexp = 1;
endexp = 20001;
classifierparamsfiles = {'../JAABA_classifier_params1.txt','../JAABA_classifier_params2.txt','../JAABA_classifier_params3.txt'};

classifierfiles = {};
configfiles = {};
for k = 1:numel(classifierparamsfiles),
classifierparamsfile = classifierparamsfiles{k};

fid = fopen(classifierparamsfile,'r');
while true,
  l = fgetl(fid);
  if ~ischar(l),
    break;
  end
  if isempty(l),
    continue;
  end
  ws = regexp(l,',','split');
  classifierfiles{end+1} = ws{1};  %#ok<SAGROW>
  configfiles{end+1} = ws{2};  %#ok<SAGROW>
end
fclose(fid);
end

scorefiles = cell(1,numel(configfiles));
for i = 1:numel(configfiles),
  
  % read parameters
  fprintf('Reading classifier parameters...\n');
  configparams = ReadXMLParams(configfiles{i});
  
  
  if ~isfield(configparams.file,'scorefilename'),
    if ~isfield(configparams,'behaviors'),
      error('configparams %d does not have field behaviors',i);
    elseif ~isfield(configparams.behaviors,'names'),
      error('configparams{%d}.behaviors does not have field names',i);
    else
      scorefilename = ['scores_',configparams.behaviors.names,'.mat'];
    end
  else
    scorefilename = configparams.file.scorefilename;
  end
  
  scorefiles{i} = scorefilename;
  
end

isok = nan(1,numel(metadata));
[~,name] = fileparts(metadata(startexp).file_system_path);
fprintf(repmat(' ',[1,numel(name)]));
for i = 1:numel(scorefiles),
  tmp = regexpi(scorefiles{i},'^scores_?(.+)\.mat$','tokens','once');
  tmp = tmp{1};
  fprintf('\t%s',tmp);
  if numel(tmp) < 20,
    fprintf(repmat(' ',[1,20-numel(tmp)]));
  end
end
fprintf('\n');
for i = startexp:endexp,
  if i > numel(metadata), break; end
  isscore = false(1,numel(scorefiles));
  [~,name] = fileparts(metadata(i).file_system_path);
  fprintf('%s',name);
  for j = 1:numel(scorefiles),
    tmp = dir(fullfile(metadata(i).file_system_path,scorefiles{j}));
    isscore(j) = ~isempty(tmp);
    if ~isempty(tmp),
      fprintf('\t%s',tmp.date);
    else
      fprintf('\tMISSING');
    end
  end
  isok(i) = all(isscore);
  fprintf('\n');
end