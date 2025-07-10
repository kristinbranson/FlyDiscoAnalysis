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
    addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;
    addpath '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    %rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    rootdatadir = '/groups/branson/bransonlab/projects/olympiad/HackHitData/';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% parameters

timestamp = '20120312';
figfiledir = sprintf('fractimefigs_%s',timestamp);
analysis_protocol = 'current';
bowlhackhitmatfile = sprintf('BowlHackHitData%s.mat',timestamp);
nbins_bowl = 5;
min_nexps_per_set = 2;
bowl_fns = {'fractime_chase_male','fractime_chase_female',...
  'fractime_walk_any','fractime_walk_male','fractime_walk_female',...
  'fractime_stop_any','fractime_stop_male','fractime_stop_female',...
  'fractime_jump_any','fractime_jump_male','fractime_jump_female'};
nbowl_fns = numel(bowl_fns);
fractime_fns = {'fractime_chase_male','fractime_chase_female',...
    'fractime_walk_any','fractime_walk_male','fractime_walk_female',...
  'fractime_stop_any','fractime_stop_male','fractime_stop_female',...
  'fractime_jump_any','fractime_jump_male','fractime_jump_female'};
nfractime_fns = numel(fractime_fns);
compute_hit_probability_params = {'nbins_control',5};
%hackhitsparamsfilestr = 'hack_hit_category_params.txt';
%hackhitsparamsfile = fullfile(settingsdir,analysis_protocol,hackhitsparamsfilestr);
%params = ReadParams(hackhitsparamsfile);

% alice's thresholds
bowlthresh = struct;
bowlthresh.fractime_moving_low = (.39+.42)/2;
bowlthresh.fractime_moving_high = (.67+.73)/2;
bowlthresh.meanspeed_moving_low = (10+11.5)/2;
bowlthresh.meanspeed_moving_high = (14.25+15)/2;
%bowlthresh.fractime_jumping_low = (9e-05+2e-04)/2;
%bowlthresh.fractime_jumping_high = ( 0.0011+0.0014)/2;
bowlthresh.fractime_close_low = (.08+.11)/2;
bowlthresh.fractime_close_high = (.23+.4)/2;
% bowlthresh.fractime_chase_female_low = 2.0726e-04;
% bowlthresh.fractime_chase_female_high = .004;
% bowlthresh.fractime_chase_male_low = 9.8378e-04;
% bowlthresh.fractime_chase_male_high = .0909;

DOPLOTRAWDATA = true;
DOPLOTSETS = true;
DOPLOTPHIT = false;
DOPLOTSORTEDDATA = true;
DOHISTOGRAMDATA = true;

%SAGEParams = {'checkflags',true,'removemissingdata',true,'screen_type','non_olympiad_mushroombody'};
%SAGEParams = {'checkflags',true,'removemissingdata',true,'screen_type','primary'};
SAGEParams = {'screen_type','primary','daterange',{'20120222T000000','20120224T000000'},...
  'checkflags',true,'removemissingdata',true};
% line_names = {'GMR_31B05_AE_01'
%   'GMR_26B07_AE_01'
%   'GMR_23B09_AE_01'
%   'GMR_47D11_AE_01'
%   'GMR_64G05_AE_01'
%   'GMR_48G07_AE_01'
%   'GMR_72H06_AE_01'
%   'GMR_81B10_AE_01'
%   'GMR_88C08_AE_01'
%   'GMR_31A11_AE_01'
%   'GMR_26B10_AE_01'
%   'GMR_26A06_AE_01'
%   'pBDPGAL4U'};
% SAGEParams = {'screen_type','primary','line_name',line_names,'checkflags',true,'removemissingdata',true};

%% load in all the FlyBowl data

if ~exist(bowlhackhitmatfile,'file'),
  metadata = SAGEListBowlExperiments(SAGEParams{:});
  metadata_exps = {metadata.experiment_name};
  metadata_exps = cellfun(@(s) s(numel('FlyBowl_')+1:end),metadata_exps,'UniformOutput',false);
  bowldata = LoadHackHitFlyBowlData('settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol,'rootdatadir',rootdatadir,...
    'experiment_names',metadata_exps);
  [~,idx] = ismember({bowldata.experiment_name},metadata_exps);
  fns_metadata = fieldnames(metadata);
  fns_bowldata = fieldnames(bowldata);
  fnsadd = setdiff(fns_metadata,fns_bowldata);
  for i = 1:numel(bowldata),
    for j = 1:numel(fnsadd),
      bowldata(i).(fnsadd{j}) = metadata(idx(i)).(fnsadd{j});
    end
  end
  save(bowlhackhitmatfile,'bowldata');
else
  load(bowlhackhitmatfile,'bowldata');
end

nbowldata = numel(bowldata);

% check for completeness
for i = 1:nbowl_fns,
  fn = ['mean_',bowl_fns{i}];
  if numel([bowldata.(fn)]) ~= nbowldata,
    error('Number of values for %s does not match number of entries in bowldata',fn);
  end
end


%% compute hit probabilities

data = struct;
data.line_names = {bowldata.line_name};
data.experiment_names = {bowldata.set};
data.session_names = {bowldata.experiment_name};
bowlphit = struct;
for i = 1:nbowl_fns,
  fn = ['mean_',bowl_fns{i}];
  data.x = [bowldata.(fn)];
  for dothreshhigh = [false,true],
    if dothreshhigh,
      fnout = [bowl_fns{i},'_high'];
    else
      fnout = [bowl_fns{i},'_low'];
    end
    if isfield(bowlthresh,fnout),
      params = {'session_thresh',bowlthresh.(fnout)};
    else
      params = {'session_falsepositiverate',.025};
    end
    bowlphit.(fnout) = ComputeHitProbability(data,'dothreshhigh',dothreshhigh,...
      params{:},compute_hit_probability_params{:});
  end
end

%% add raw statistics to bowlphit

[line_names,~,lineidx] = unique({bowldata.line_name});
nlines = numel(line_names);

% compute per-line statistics
bowl_linedata = struct;

for i = 1:nlines,
  idx = lineidx == i;
  bowl_linedata_curr = struct;
  bowl_linedata_curr.line_name = line_names{i};
  for j = 1:nbowl_fns,
    fn = bowl_fns{j};
    bowl_linedata_curr.(fn) = struct;
    tmp = [bowldata(idx).(['mean_',fn])];
    bowl_linedata_curr.(fn).mean = nanmean(tmp);
    bowl_linedata_curr.(fn).std = nanstd(tmp,1);
    bowl_linedata_curr.(fn).n = nnz(~isnan(tmp));
  end
  if i == 1,
    bowl_linedata = repmat(bowl_linedata_curr,[1,nlines]);
  else
    bowl_linedata(i) = bowl_linedata_curr;
  end
end

% output to file
outfilename = sprintf('BowlLineData_%s.csv',timestamp);
fid = fopen(outfilename,'w');
fprintf(fid,'line_name');
for j = 1:nbowl_fns,
  fn = bowl_fns{j};
  fprintf(fid,',%s_mean,%s_std,%s_n',fn,fn,fn);
end
fprintf(fid,'\n');

for i = 1:nlines,
  s = [];
  fprintf(fid,'%s',bowl_linedata(i).line_name);
  for j = 1:nbowl_fns,
    fn = bowl_fns{j};
    fprintf(fid,',%f,%f,%f',bowl_linedata(i).(fn).mean,...
      bowl_linedata(i).(fn).std,bowl_linedata(i).(fn).n);
  end
  fprintf(fid,'\n');
end

fclose(fid);

for j = 1:nbowl_fns,
  fn = bowl_fns{j};
  high_fn = [fn,'_high'];
  low_fn = [fn,'_high'];
  [ism,tmplineidx] = ismember(line_names,bowlphit.(high_fn).line.names);
  bowlphit.(high_fn).line.mean = nan(size(bowlphit.(high_fn).line.names));
  bowlphit.(high_fn).line.std = nan(size(bowlphit.(high_fn).line.names));
  for i = find(ism),
    bowlphit.(high_fn).line.mean(tmplineidx(i)) = bowl_linedata(i).(fn).mean;
    bowlphit.(high_fn).line.std(tmplineidx(i)) = bowl_linedata(i).(fn).std;
  end
  [ism,tmplineidx] = ismember(line_names,bowlphit.(low_fn).line.names);
  bowlphit.(low_fn).line.mean = nan(size(bowlphit.(low_fn).line.names));
  bowlphit.(low_fn).line.std = nan(size(bowlphit.(low_fn).line.names));
  for i = find(ism),
    bowlphit.(low_fn).line.mean(tmplineidx(i)) = bowl_linedata(i).(fn).mean;
    bowlphit.(low_fn).line.std(tmplineidx(i)) = bowl_linedata(i).(fn).std;
  end
end


% compute per-set statistics
[set_names,~,setidx] = unique({bowldata.set});
nsets = numel(set_names);
bowl_setdata = struct;

for i = 1:nsets,
  idx = setidx == i;
  bowl_setdata_curr = struct;
  bowl_setdata_curr.set = set_names{i};
  for j = 1:nbowl_fns,
    fn = bowl_fns{j};
    bowl_setdata_curr.(fn) = struct;
    tmp = [bowldata(idx).(['mean_',fn])];
    bowl_setdata_curr.(fn).mean = nanmean(tmp);
    bowl_setdata_curr.(fn).std = nanstd(tmp,1);
    bowl_setdata_curr.(fn).n = nnz(~isnan(tmp));
  end
  if i == 1,
    bowl_setdata = repmat(bowl_setdata_curr,[1,nsets]);
  else
    bowl_setdata(i) = bowl_setdata_curr;
  end
end


for j = 1:nbowl_fns,
  fn = bowl_fns{j};
  high_fn = [fn,'_high'];
  low_fn = [fn,'_low'];
  [ism,tmpsetidx] = ismember(set_names,bowlphit.(high_fn).exp.names);
  bowlphit.(high_fn).exp.mean = nan(size(bowlphit.(high_fn).exp.names));
  bowlphit.(high_fn).exp.std = nan(size(bowlphit.(high_fn).exp.names));
  for i = find(ism),
    bowlphit.(high_fn).exp.mean(tmpsetidx(i)) = bowl_setdata(i).(fn).mean;
    bowlphit.(high_fn).exp.std(tmpsetidx(i)) = bowl_setdata(i).(fn).std;
  end
  [ism,tmpsetidx] = ismember(set_names,bowlphit.(low_fn).exp.names);
  bowlphit.(low_fn).exp.mean = nan(size(bowlphit.(low_fn).exp.names));
  bowlphit.(low_fn).exp.std = nan(size(bowlphit.(low_fn).exp.names));
  for i = find(ism),
    bowlphit.(low_fn).exp.mean(tmpsetidx(i)) = bowl_setdata(i).(fn).mean;
    bowlphit.(low_fn).exp.std(tmpsetidx(i)) = bowl_setdata(i).(fn).std;
  end
end

%% save results to file

save(bowlhackhitmatfile,'bowlphit','bowldata');

%% combine data per line

[line_names,~,lineidx] = unique({bowldata.line_name});
nlines = numel(line_names);

% compute the means and standard deviations for each line
nexpsperline = zeros(1,nlines);
for j = 1:nlines,
  idx = lineidx == j;
  if ~any(idx),
    continue;
  end
  nexpsperline(j) = nnz(idx);
end
[nexpsperline,order] = sort(nexpsperline); %#ok<ASGLU>
line_names = line_names(order);
[~,lineidx] = ismember({bowldata.line},line_names);

line_names_print = regexprep(line_names,'^GMR_','');
line_names_print = regexprep(line_names_print,'_01$','');
line_names_print = regexprep(line_names_print,'_AE$','');
line_names_print = regexprep(line_names_print,'_AD$','D');


controlidx = find(strcmpi(line_names,'pBDPGAL4U'));
bowlmu = nan(nbowl_fns,nlines);
bowlsig = nan(nbowl_fns,nlines);
nexpsperline = zeros(1,nlines);
for i = 1:nbowl_fns,
  fn = ['mean_',bowl_fns{i}];
  for j = 1:nlines,
    idx = lineidx == j;
    if ~any(idx),
      continue;
    end
    bowlmu(i,j) = nanmean([bowldata(idx).(fn)]);
    bowlsig(i,j) = nanstd([bowldata(idx).(fn)],1);
    if i == 1,
      nexpsperline(j) = nnz(idx);
    end
  end
end

if ~isempty(figfiledir) && ~exist(figfiledir,'dir'),
  mkdir(figfiledir);
end

%% plot the raw data

if DOPLOTRAWDATA,
  
% plot the bowl data
hfig = 1;
figure(hfig);
clf;
%hax = nan(1,nbowl_fns);
colors = jet(10)*.7;
hax = createsubplots(nbowl_fns,1,[[.025,.025];[.05,.025]]);
for i = 1:nbowl_fns,
  fn = ['mean_',bowl_fns{i}];
  axes(hax(i)); %#ok<LAXES>
  plot(lineidx,[bowldata.(fn)],'k.');
  hold on;
  for j = 1:size(colors,1),
    idx = (j-1)+1:size(colors,1):nlines;
    plot([idx;idx],[bowlmu(i,idx)-bowlsig(i,idx);bowlmu(i,idx)+bowlsig(i,idx)],'-','color',colors(j,:));
    plot(idx,bowlmu(i,idx),'.','color',colors(j,:));
  end
  hti = title(fn);
  set(hti,'interpreter','none');
  axisalmosttight;
  set(hax(i),'XTickLabel',{});
  plot([.5,nlines+.5],bowlmu(i,controlidx)+[0,0],'c--');
  hold on;
  plot([.5,nlines+.5],bowlmu(i,controlidx)+bowlsig(i,controlidx)+[0,0],'c--');
  plot([.5,nlines+.5],bowlmu(i,controlidx)-bowlsig(i,controlidx)+[0,0],'c--');  
end
set(hax,'XLim',[.5,nlines+.5]);
linkaxes(hax,'x');
xtick = unique(round(linspace(1,nlines,50)));
set(hax,'XTick',xtick);
set(hax(end),'XTickLabel',line_names_print(xtick));
th = rotateticklabel(hax(end),90);  %#ok<NASGU>


set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
filename = sprintf('BowlLineVsPerSessionRawData_%s',timestamp);
SaveFigsHelper(hfig,filename,figfiledir);

end

% set(hax(1),'XLim',[nlines-49.5,nlines+.5]);
% delete(th);
% set(hax(end),'XTick',nlines-50:nlines,'XTickLabel',line_names_print(nlines-50:nlines));
% th = rotateticklabel(hax(end),90); %#ok<NASGU>
% set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
% filename = sprintf('BowlLineVsPerSessionRawData_Zoom_%s',timestamp);
% SaveFigsHelper(hfig,filename);

%% plot sets

if DOPLOTSETS,

[set_names,~,setidx] = unique({bowldata.set});
nsets = numel(set_names);
setdata = struct;
notenoughdata = false(1,nsets);
for i = 1:nsets,
  idx = setidx == i;
  if nnz(idx) < min_nexps_per_set,
    notenoughdata(i) = true;
  end
  k = find(idx,1);
  setdata(i).set_name = set_names{i};
  setdata(i).line_name = bowldata(k).line_name;
  setdata(i).idx = find(idx);
  setdata(i).nexps = nnz(idx);
  for j = 1:nbowl_fns,
    fn = ['mean_',bowl_fns{j}];
    setdata(i).(['mean',fn]) = nanmean([bowldata(idx).(fn)]);
    setdata(i).(['std',fn]) = nanstd([bowldata(idx).(fn)],1);
  end
end
setdata(notenoughdata) = [];
[set_names,~,setidx] = unique({setdata.set_name});
nsets = numel(set_names);
% setdata = ClassifyHackHitCategory(setdata,'settingsdir',settingsdir,...
%   'analysis_protocol',analysis_protocol,...
%   'issetdata',true);

[ism,setlineidx] = ismember({setdata.line_name},line_names);

hfig = 2;
figure(hfig);
clf;
bowlsetmu = nan(nbowl_fns,nlines);
bowlsetsig = nan(nbowl_fns,nlines);
nsetsperline = zeros(1,nlines);
for i = 1:nbowl_fns,
  fn = ['meanmean_',bowl_fns{i}];
  for j = 1:nlines,
    idx = setlineidx == j;
    if ~any(idx),
      continue;
    end
    bowlsetmu(i,j) = nanmean([setdata(idx).(fn)]);
    bowlsetsig(i,j) = nanstd([setdata(idx).(fn)],1);
    if i == 1,
      nsetsperline(j) = nnz(idx);
    end
  end
end

clf;
hax = createsubplots(nbowl_fns,1,[[.025,.025];[.05,.025]]);
colors = jet(10)*.7;
for i = 1:nbowl_fns,
  fn = ['meanmean_',bowl_fns{i}];
  stdfn = ['stdmean_',bowl_fns{i}];
  %hax(i) = subplot(nbowl_fns,1,i);
  axes(hax(i)); %#ok<LAXES>
  plot([setlineidx;setlineidx],[[setdata.(fn)]+[setdata.(stdfn)];[setdata.(fn)]-[setdata.(stdfn)]],'-','color',[.7,.7,.7]);
  hold on;
  plot(setlineidx,[setdata.(fn)],'k.');
  for j = 1:size(colors,1),
    idx = (j-1)+1:size(colors,1):nlines;
    plot([idx;idx],[bowlsetmu(i,idx)-bowlsetsig(i,idx);bowlsetmu(i,idx)+bowlsetsig(i,idx)],'-','color',colors(j,:));
    plot(idx,bowlsetmu(i,idx),'.','color',colors(j,:));
  end
  hti = title(fn);
  set(hti,'interpreter','none');
  axisalmosttight;
  set(hax(i),'XTickLabel',{});
  plot([1,nlines],bowlsetmu(i,controlidx)+[0,0],'c--');
  hold on;
  plot([1,nlines],bowlsetmu(i,controlidx)+bowlsetsig(i,controlidx)+[0,0],'c--');
  plot([1,nlines],bowlsetmu(i,controlidx)-bowlsetsig(i,controlidx)+[0,0],'c--');  
end
set(hax,'XLim',[1,nlines]);
linkaxes(hax,'x');
xtick = unique(round(linspace(1,nlines,50)));
set(hax,'XTick',xtick);
set(hax(end),'XTickLabel',line_names_print(xtick));
th = rotateticklabel(hax(end),90);

set(hfig,'Units','pixels','Position',[43 439 1037 1344]);
filename = sprintf('BowlLineVsPerExperimentRawData_%s',timestamp);
SaveFigsHelper(hfig,filename,figfiledir);
% 
% set(hax(1),'XLim',[nlines-49.5,nlines+.5]);
% delete(th);
% set(hax(end),'XTick',nlines-50:nlines,'XTickLabel',line_names_print(nlines-50:nlines));
% th = rotateticklabel(hax(end),90);
% set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
% filename = sprintf('BowlLineVsPerExperimentRawData_Zoom_%s',timestamp);
% SaveFigsHelper(hfig,filename);

end

%% histogram the raw data

if DOHISTOGRAMDATA,

% choose bins for each stat
bowlminv = nan(1,nbowl_fns);
bowlmaxv = nan(1,nbowl_fns);
for i = 1:nbowl_fns,
  fn = ['meanmean_',bowl_fns{i}];
  bowlminv(i) = min([setdata.(fn)]);
  bowlmaxv(i) = max([setdata.(fn)]);
end

hfig = 3;
figure(hfig);
clf;
w = zeros(1,nsets);
for i = 1:nlines,
  idx = setlineidx == i;
  if ~any(idx),
    continue;
  end
  w(idx) = 1/nsetsperline(i);
end
hax = createsubplots(nbowl_fns,1,[[.025,.025];[.04,.04]]);

for i = 1:nbowl_fns,
  fn = ['meanmean_',bowl_fns{i}];
  axes(hax(i)); %#ok<LAXES>
  [counts,ctrs] = myhist([setdata.(fn)],20,'weights',w);
  plot(ctrs,counts/nlines,'k.-');
  hold on;
  ylim = get(hax(i),'YLim');
  plot(bowlsetmu(i,controlidx)+[0,0],ylim,'c--');
  plot(bowlsetmu(i,controlidx)+bowlsetsig(i,controlidx)+[0,0],ylim,'c--');
  plot(bowlsetmu(i,controlidx)-bowlsetsig(i,controlidx)+[0,0],ylim,'c--');  
  hti = title(bowl_fns{i});
  set(hti,'Interpreter','none');
end

set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
filename = sprintf('BowlSessionRawDataHistogram_%s',timestamp);
SaveFigsHelper(hfig,filename,figfiledir);

end

%% plot sorted set data 

if DOPLOTSORTEDDATA,

hfig = 4;
figure(hfig);
clf;
hax = createsubplots(nbowl_fns,1,[[.025,.025];[.05,.025]]); 
colors = jet(10)*.7;
for i = 1:nbowl_fns,

  isbaddata = isnan(bowlsetmu(i,:));
  [~,order] = sort(bowlsetmu(i,~isbaddata));
  gooddataidx = find(~isbaddata);
  baddataidx = find(isbaddata);
  order = [baddataidx,gooddataidx(order)]; 
  sorted_bowlsetmu = bowlsetmu(i,order);
  sorted_bowlsetsig = bowlsetsig(i,order);
  sorted_linenames = line_names(order);
  [~,sorted_setlineidx] = ismember({setdata.line_name},sorted_linenames);
  
  fn = ['meanmean_',bowl_fns{i}];
  stdfn = ['stdmean_',bowl_fns{i}];
  axes(hax(i)); %#ok<LAXES>
  plot([sorted_setlineidx;sorted_setlineidx],[[setdata.(fn)]+[setdata.(stdfn)];[setdata.(fn)]-[setdata.(stdfn)]],'-','color',[.7,.7,.7]);
  hold on;
  plot(sorted_setlineidx,[setdata.(fn)],'k.');
  for j = 1:size(colors,1),
    idx = (j-1)+1:size(colors,1):nlines;
    plot([idx;idx],[sorted_bowlsetmu(idx)-sorted_bowlsetsig(idx);sorted_bowlsetmu(idx)+sorted_bowlsetsig(idx)],'-','color',colors(j,:));
    plot(idx,sorted_bowlsetmu(idx),'.','color',colors(j,:));
  end
  hti = title(bowl_fns{i});
  set(hti,'interpreter','none');
  axisalmosttight;
  sorted_controlidx = find(strcmpi(sorted_linenames,'pBDPGAL4U'),1);
  set(hax(i),'XTick',sorted_controlidx,'XTickLabel','pBDPGAL4U');
  plot([1,nlines],bowlsetmu(i,controlidx)+[0,0],'c--');
  hold on;
  plot([1,nlines],bowlsetmu(i,controlidx)+bowlsetsig(i,controlidx)+[0,0],'c--');
  plot([1,nlines],bowlsetmu(i,controlidx)-bowlsetsig(i,controlidx)+[0,0],'c--');  
end
set(hax,'XLim',[1,nlines]);


set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
filename = sprintf('BowlLineVsPerExperimentRawData_Sorted_%s',timestamp);
SaveFigsHelper(hfig,filename,figfiledir);
% 
% set(hax,'XLim',[nlines-49.5,nlines]);
% set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
% filename = sprintf('BowlLineVsPerExperimentRawData_Sorted_Zoom_%s',timestamp);
% SaveFigsHelper(hfig,filename);

end

%% plot phit

if DOPLOTPHIT,

bowlphit_fns = fieldnames(bowlphit); %#ok<*UNRCH>
nbowlphit_fns = numel(bowlphit_fns);
hax = nan(1,nbowlphit_fns);
colors = jet(10)*.7;

% sort lines by number of retests
line_names = {};
for i = 1:nbowlphit_fns,
  line_names = union(line_names,bowlphit.(bowlphit_fns{i}).line.names);
end
nlines = numel(line_names);
nexpsperline = zeros(1,nlines);
for i = 1:nbowlphit_fns,
  [~,idx] = ismember(bowlphit.(bowlphit_fns{i}).line.names,line_names);
  nexpsperline(idx) = nexpsperline(idx) + bowlphit.(bowlphit_fns{i}).line.n;
end
[~,lineorder] = sort(nexpsperline);
line_names = line_names(lineorder);

%hax = createsubplots(ceil(nbowlphit_fns/2),2,.05);
%hax = reshape(hax,[ceil(nbowlphit_fns/2),2])';

epsilon = .5;
for i = 1:nbowlphit_fns,
  hfig = 6+i;
  figure(hfig);
  clf;
 
  %hax(i) = subplot(ceil(nbowlphit_fns/2),2,i);
  fn = bowlphit_fns{i};

%   % plot per-session data
%   [~,session_lineidx] = ismember(bowlphit.(fn).session.line_names,line_names);
%   minx = min(bowlphit.(fn).session.x);
%   maxx = max(bowlphit.(fn).session.x);
%   plot(session_lineidx,(bowlphit.(fn).session.x-minx)/(maxx-minx),'.','color',[.7,.7,.7]);
%  hold on;

  hax(i) = gca;
  hold on;

  exp_jitter = epsilon*(rand(1,numel(bowlphit.(fn).exp.frachits))-.5);
  % plot per-exp data
  [~,exp_lineidx] = ismember(bowlphit.(fn).exp.line_names,line_names);
  plot(exp_lineidx+exp_jitter,bowlphit.(fn).exp.frachits,'.','color',[.7,.7,.7]);

  % plot per-line data
  [~,line_lineidx] = ismember(bowlphit.(fn).line.names,line_names);
  for j = 1:size(colors,1),
    idx = mod(line_lineidx,size(colors,1)) == j-1;
    plot(line_lineidx(idx),bowlphit.(fn).line.phit(idx),'o','color',colors(j,:),'markerfacecolor',colors(j,:));
  end
  
  % plot per-exp data
  for j = 1:size(colors,1),
    idx = mod(exp_lineidx,size(colors,1)) == j-1;
    plot(exp_lineidx(idx)+exp_jitter(idx),bowlphit.(fn).exp.phit(idx),'.','color',(colors(j,:)+.7)/2);
  end
  
  hti = title(fn);
  set(hti,'interpreter','none');
  axisalmosttight;
  %set(hax(i),'XTickLabel',{});
end
set(hax,'XLim',[1,nlines]);
linkaxes(hax,'x');
xtick = unique(round(linspace(1,nlines,50)));
set(hax,'XTick',xtick);
set(hax,'XTickLabel',line_names_print(xtick));
th = cell(1,nbowlphit_fns);
for i = 1:nbowlphit_fns,
  th{i} = rotateticklabel(hax(i),90);
end
drawnow;

for i = 1:nbowlphit_fns,
  hfig = 6+i;
  set(hax(i),'XGrid','on','XLim',[.5,nlines+.5]);
  set(hfig,'Units','pixels','Position',[43 439 1037 1003]);
  filename = sprintf('BowlLineVsPHit_%s_%s',bowlphit_fns{i},timestamp);
  SaveFigsHelper(hfig,filename,figfiledir);
  
end

end