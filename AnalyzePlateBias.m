%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;

%% parameters

rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
analysis_protocol = '20110211';
datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
biasdiagnosticsparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.biasdiagnosticsparamsfilestr);
params = ReadParams(biasdiagnosticsparamsfile);
registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
registration_params = ReadParams(registrationparamsfile);

%% get experiment directories

data = SAGEGetBowlData('daterange',{'20110216T000000','20110217T999999'},...
  'plate',{'11','12','13'},'data_type','bias_diagnostics*',...
  'experiment_name','FlyBowl_*');

%% get data for plates 10, 14

data_flat = SAGEGetBowlData('daterange',{'20110216T000000','20110217T999999'},...
  'plate',{'10','14'},'data_type','bias_diagnostics*',...
  'experiment_name','FlyBowl_*');

%% show directories

experiment_names = {data.experiment_name};
experiment_names = cellfun(@(s) s(9:end),experiment_names,'UniformOutput',false);
experiment_names = sort(experiment_names);
for i = 1:numel(experiment_names),
  web(fullfile(rootdatadir,experiment_names{i}),'-browser');
end

%% combine all data

dataall = structappend(data,data_flat);

%% compute maximum radius

[edges_r,centers_r] = SelectHistEdges(params.nbins_r,[0,1]*registration_params. circleRadius_mm,'linear');
nexps = numel(dataall);
for i = 1:nexps,
  fracsmooth = dataall(i).bias_diagnostics_fracsmooth;
  fracsmooth = reshape(fracsmooth,[params.nbins_r,params.nbins_theta]);
  [maxfracsmooth,maxr_bin] = max(fracsmooth,[],1);
  dataall(i).argmaxr_fracsmooth = centers_r(maxr_bin);
  dataall(i).r_farthest = max(dataall(i).argmaxr_fracsmooth);
end

%% plot r_farthest

plate = [dataall.plate];
bowl = {dataall.bowl};
[allplates,~,idx_plate] = unique(plate);
[allbowls,~,idx_bowl] = unique(bowl);
idx = (idx_plate-1)*numel(allbowls)+(idx_bowl-1);

figure(1);
clf;
hax = gca;
plot(idx,[dataall.r_farthest],'.');
s = cell(numel(allbowls),numel(allplates));
for platei = 1:numel(allplates),
  for bowli = 1:numel(allbowls),
    s{bowli,platei} = sprintf('%d%s',allplates(platei),allbowls{bowli});
  end
end
set(hax,'XTick',0:numel(allplates)*numel(allbowls)-1,'XTickLabel',s(:));