% ScriptCheckExperiments20120726

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/starvation;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/age;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211/results';
analysis_protocol = 'current_non_olympiad_azanchir_mating_galit_CS_20120211';
mintimestamp = '20120501T000000';
npar = 8;


%% set data locations

metadata = SAGEListBowlExperiments('checkflags',false,'removemissingdata',false,...
  'handling_protocol',{'HP_flybowl_v007p1.xls','HP_flybowl_v007p2.xls','HP_flybowl_v007p3.xls',...
  'HP_flybowl_v007p6.xls','HP_flybowl_v007p7.xls'},...
  'screen_type','non_olympiad_azanchir_mating_galit_CS_20120211',...
  'rootdir',rootdatadir);
% fid = fopen(fullfile(outdatadir,outexpdirfilename),'w');
% fprintf(fid,'%s\n',metadata.file_system_path);
% fclose(fid);


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
    save(sprintf('CheckExperiments_galit_%s_20120910.mat',datatype),'issuccess','msgs','order','ii','metadata');
  end
  
end

%% check for scores files

startexp = 1;
endexp = 20001;
classifierparamsfiles = {'../JAABA_classifier_params1.txt','../JAABA_classifier_params2.txt'};

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