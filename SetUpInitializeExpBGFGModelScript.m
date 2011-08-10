% set up for LearnExpBgFGModelScript.py

%% set up path, data locations

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
if ispc,  
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  ctrax_settingsdir = 'E:\Code\FlyBowlCtrax';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  ctrax_settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax';
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/00_incoming';
  rootdir_fixed = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/Fixed_AR/EP00005_rc8';
  %rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
end
analysis_protocol = '20110804';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};
annfilestr = 'testmovie.ufmf.ann';

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

tmpdir = ctrax_settingsdir;
if ~exist(tmpdir,'file'),
  mkdir(tmpdir);
end
tmpdir = fullfile(ctrax_settingsdir,analysis_protocol);
if ~exist(tmpdir,'file'),
  mkdir(tmpdir);
end
paramsdir = tmpdir;
tmpdir = fullfile(tmpdir,'LearnCtraxParams');
if ~exist(tmpdir,'file'),
  mkdir(tmpdir);
end
outputdir = tmpdir;
expdirsFileStr = 'expdirs.txt';
expdirsTestFileStr = 'expdirs_test.txt';

%% experiments to analyze

experiment_params = struct;
% experiment_params.experiment_name = 'FlyBowl_*';
% % check for failures
% experiment_params.checkflags = false;
% % rearing condition
% experiment_params.rearing_protocol = 'RP_Olympiad_v008p2.xls';
% % experimental protocol
% experiment_params.experiment_protocol = 'EP_flybowl_v010p0.xls';
% root data dir
experiment_params.rootdir = rootdir;
% success file
experiment_params.subreadfiles = {'SUCCESS'};

% manual_pf
%experiment_params.manual_pf = {'P','U'};
% data types
%experiment_params.data_type = {'registrationdata_circleCenterX',...
%  'registrationdata_circleCenterY',...
%  'registrationdata_circleRadius'};
%experiment_params.dataset = 'score';
ngal4 = 0;
ncontrols = 8;
tmpparams = struct2paramscell(experiment_params);

% output file
expdirsfilename = fullfile(outputdir,expdirsFileStr);
expdirstestfilename = fullfile(outputdir,expdirsTestFileStr);

%experiments_all = SAGEListBowlExperiments(tmpparams{:});
[expdirs_all,expdir_reads_all,expdir_writes_all,experiments_all] = getExperimentDirs(tmpparams{:});
%[experiments_all.line] = deal(experiments_all.line_name);

%% choose some experiments to track

[experiments,ngal4_chosen,ncontrols_chosen] = choose_expdirs(experiments_all,ngal4,ncontrols);

%% compute background models for all experiments

fprintf('Run this command from the Ctrax directory to set up initial tracking parameters:\n');
fprintf('python Ctrax.py --Input=%s/movie.ufmf --Output=%s/%s --SettingsFile=%s/settings.ann\n',...
  experiments(1).file_system_path,experiments(1).file_system_path,annfilestr,fullfile(ctrax_settingsdir,'current'));
fprintf('Let it track for a while to generate a settings file for use in other movies.\n');
settingsfile = fullfile(experiments(1).file_system_path,'testmovie.ufmf.ann');

expdirfile = fullfile(ctrax_settingsdir,analysis_protocol,'LearnCtraxParams','expdirs_initialize.txt');
fid = fopen(expdirfile,'w');
if fid < 0,
  error('Could not open file %s for writing',expdirfile);
end
fprintf(fid,'%s\n',experiments.file_system_path);
fclose(fid);
mkdir(fullfile(ctrax_settingsdir,analysis_protocol,'LearnCtraxParams'),'out');

fprintf('After tracking is done, run tracking for maybe 20 minutes with settings file %s on all of the experiment directories listed in %s:\n',settingsfile,expdirfile);

% ./runall.pl /groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/expdirs_initialize.txt ./runctrax20110804.sh /groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/out

%% look at background models for all experiments

isgoodbg = false(1,numel(experiments));

for iscontrol = [false,true],

  idx = find(strcmpi({experiments.line},'pBDPGAL4U') == iscontrol);
  if iscontrol,
    maxn = ncontrols;
    fprintf('Choosing controls...\n');
  else
    maxn = ngal4;
    fprintf('Choosing GAL4s...\n');
  end
  n = 0;
  for ii = 1:numel(idx),
    i = idx(ii);
   
    exp = experiments(i);
    annfilename = fullfile(exp.file_system_path,annfilestr);
    moviefilename = fullfile(exp.file_system_path,dataloc_params.moviefilestr);
    [bgcenter,bgstd,movie_height,movie_width] = ...
      read_ann(annfilename,'background_center','background_mad','movie_height','movie_width');
    if isempty(bgcenter),
      error('Error reading ann file %s',annfilename);
    end
    if isempty(movie_height) || isempty(movie_width),
      [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
      movie_height = headerinfo.nr;
      movie_width = headerinfo.nc;
      fclose(fid);
    end

    bgcenter = permute(reshape(bgcenter,[movie_height,movie_width,numel(bgcenter)/...
      (movie_height*movie_width)]),[2,1,3]);
    bgstd = permute(reshape(bgstd,[movie_height,movie_width,numel(bgstd)/...
      (movie_height*movie_width)]),[2,1,3]);
    hold off; 
    subplot(1,2,1);
    imagesc(bgcenter,[0,255]); 
    axis image;
    subplot(1,2,2);
    imagesc(bgstd); axis image;
    
    isgoodbg(i) = input('Enter 1 if good, 0 if bad: ');
    
    if isgoodbg(i),
      n = n + 1;
    end
    if n >= maxn,
      break;
    end
    
  end
    
end

%% track all these experiments

fprintf('Run tracking with settings file %s on all of the experiment directories listed in %s:\n',settingsfile,expdirfile);

% ./runall.pl /groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/expdirs_initialize.txt ./runctrax20110804.sh /groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/out

%% output chosen experiments to file

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);
fid = fopen(expdirsfilename,'w');
for i = 1:nexpdirs,
  fprintf(fid,'%s\n',expdirs{i});
end
fclose(fid);

%% command for learning

paramsFileStr = 'ExpBGFGModelParams.txt';
movieFileStr = dataloc_params.moviefilestr;
annFileStr = annfilestr;
outputFileStr = 'ExpBGFGModelResults.pickle';
matFileStr = 'ExpBGFGModelResults.mat';

% file to write experiment directory names to
expdirsFileName = fullfile(outputdir,expdirsFileStr);
expdirsFileNameTest = fullfile(outputdir,expdirsTestFileStr);

% file containing parameters
paramsFileName = fullfile(paramsdir,paramsFileStr);

% file to write results to
outputFileName = fullfile(outputdir,outputFileStr);

% mat file to write results to
matFileName = fullfile(outputdir,matFileStr);

fprintf('execute the following command:\n');
fprintf('python ExpBGFGModel.py');
fprintf(' -f %s',expdirsFileName);
fprintf(' -p %s',paramsFileName);
fprintf(' -m %s',movieFileStr);
fprintf(' -a %s',annFileStr);
fprintf(' -o %s',outputFileName);
fprintf(' --mat %s\n',matFileName);

%% choose LoG filter threshold

fil = fspecial('log',21,5);
lim = [-5,5];
nbins = 200;
nframessample = 20;
[edges,centers] = SelectHistEdges(nbins,lim,'linear');

meanfrac = zeros(1,nbins);
foredata = [];

for expdiri = 1:nexpdirs,
  
  expdir = expdirs{expdiri};
  annfilename = fullfile(expdir,annfilestr);
  moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
  trxfilename = fullfile(expdir,dataloc_params.ctraxfilestr);
  
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
  
  [bgcenter,isarena,movie_height,movie_width] = ...
    read_ann(annfilename,'background_center','isarena','movie_height','movie_width');
  if isempty(movie_height) || isempty(movie_width),
    movie_height = headerinfo.nr;
    movie_width = headerinfo.nc;
  end
  bgcenter = permute(reshape(bgcenter,[movie_height,movie_width,numel(bgcenter)/...
    (movie_height*movie_width)]),[2,1,3]);
  isarena = permute(reshape(isarena,[movie_height,movie_width,numel(isarena)/...
    (movie_height*movie_width)]),[2,1,3]);
  res = imfilter(bgcenter,fil);
  x = res(isarena==0);
  counts = hist(x,centers);
  frac = counts / sum(counts);
  meanfrac = meanfrac + frac;
  
  trx = load_tracks(trxfilename);
  firstframe = min([trx.firstframe]);
  endframe = max([trx.endframe]);
  framessample = unique(round(linspace(firstframe,endframe,nframessample)));
  for i = 1:numel(framessample),
    f = framessample(i);
    im = readframe(f);
    res1 = imfilter(double(im),fil);
    for j = 1:numel(trx),
      if f < trx(j).firstframe || f > trx(j).endframe,
        continue;
      end
      k = f+trx(j).off;
      params = [trx(j).x(k),trx(j).y(k),trx(j).a(k)*4,trx(j).b(k)*4,trx(j).theta(k)];
      bb = [1,movie_height,1,movie_width];
      bw = ellipsepixels(params,bb,true);
      maxv = max(res1(bw));
      foredata(end+1) = maxv; %#ok<SAGROW>
    end
    
  end
  
  hold off;
  plot(centers,meanfrac/expdiri,'k.-');
  hold on;
  counts = hist(foredata,centers);
  plot(centers,counts/sum(counts),'r.-');
  title(num2str(expdiri));
  drawnow;
  
  fclose(fid);
  
end

meanfrac = meanfrac / nexpdirs;
meanfrac_fore = hist(foredata,centers);
meanfrac_fore = meanfrac_fore / sum(meanfrac_fore);

hold off;
score = cumsum(meanfrac-meanfrac_fore);
plot(centers,score,'k.-');

maxscore = max(score);
is = find(score==maxscore);
i = floor(median(is));
LoG_thresh = edges(i+1);
fprintf('LoG_thresh = %f\n',LoG_thresh);


% LoG_thresh = 1.700000

%% choose llr thresholds

n_bg_std_thresh_low = 20;
n_bg_std_thresh_high = 100;

expbgfgmodelfile = fullfile(settingsdir,analysis_protocol,dataloc_params.expbgfgmodelmatfilestr);
model = load(expbgfgmodelfile);

llr_back_store = [];
llr_fore_large_store = [];
llr_fore_small_store = [];
llr_back_edge_store = [];
llr_fore_prctiles_store = [];

lim = [-40,120];
nbins = 200;
nframessample = 20;
[edges,centers] = SelectHistEdges(nbins,lim,'linear');
[x,y] = meshgrid(1:movie_width,1:movie_height);
isnearedge = sqrt((x-movie_width/2).^2 + (y-movie_height/2).^2) >= .47*movie_width;
prctiles = [10,25,40];
colors = jet(5+numel(prctiles))*.7;
prctilestrs = cellstr(num2str(prctiles(:)));

for expdiri = 1:nexpdirs,
  
  expdir = expdirs{expdiri};
  annfilename = fullfile(expdir,annfilestr);
  moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
  trxfilename = fullfile(expdir,dataloc_params.ctraxfilestr);
  
  fprintf('%d: %s\n',expdiri,expdir);
  
  ann = read_ann(annfilename);
  ann.n_bg_std_thresh_low = n_bg_std_thresh_low;
  ann.n_bg_std_thresh_high = n_bg_std_thresh_high;
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);

  % resize images
  annfile_images = {'background_median','background_mean',...
    'background_mad','background_std',...
    'fracframesisback',...
    'background_center','background_dev',...
    'isarena'};
  
  if ~all(isfield(ann,{'movie_width','movie_height'})),
    nr = headerinfo.nr;
    nc = headerinfo.nc;
  else
    nr = ann.movie_width;
    nc = ann.movie_height;
  end
  
  % resize images read from annotation
  for i = 1:numel(annfile_images),
    fn = annfile_images{i};
    if isfield(ann,fn),
      ann.(fn) = permute(reshape(ann.(fn),[nc,nr,numel(ann.(fn))/(nr*nc)]),[2,1,3]);
    end
  end
  
  llr = ExpBGFGModel_lik(model,ann.background_center);
  llr(model.always_bg_mask == 1) = nan;
  
  thresh = prctile(llr(:),90);
  llr_back_store = [llr_back_store,llr(llr >= thresh)]; %#ok<AGROW>
  llr_back_edge_store = [llr_back_edge_store,llr(isnearedge)]; %#ok<AGROW>
  
  trx = load_tracks(trxfilename);
  firstframe = min([trx.firstframe]);
  endframe = max([trx.endframe]);
  framessample = unique(round(linspace(firstframe,endframe,nframessample)));
  for i = 1:numel(framessample),
    t = framessample(i);
    im = double(readframe(t));

    % background subtraction
    [isfore,diffim] = BackSub(im,ann);
    % use GMM with the tracked fly positions as initialization
    [cc,nfliescurr] = AssignPixels(isfore,diffim,trx,t);
    
    if ~exist('model','var'),
      llr = zeros(nr,nc);
    else
      llr = ExpBGFGModel_lik(model,im);
    end

    for k = 1:nfliescurr,
      llrcurr = llr(cc==k);
      llr_fore_large_store(end+1) = max(llrcurr); %#ok<SAGROW>
      llr_fore_small_store(end+1) = min(llrcurr); %#ok<SAGROW>
      llr_fore_prctiles_store(:,end+1) = prctile(llrcurr,prctiles); %#ok<SAGROW>
    end


    
  end
  
  hold off;
  counts = hist(llr_back_store,centers);
  plot(centers,counts/sum(counts),'.-','color',colors(1,:));
  hold on;
  counts = hist(llr_back_edge_store,centers);
  plot(centers,counts/sum(counts),'.-','color',colors(2,:));
  counts = hist(llr_fore_large_store,centers);
  plot(centers,counts/sum(counts),'.-','color',colors(3,:));
  counts = hist(llr_fore_small_store,centers);
  plot(centers,counts/sum(counts),'.-','color',colors(4,:));
  title(num2str(expdiri));
  llr_back_large_store = max(llr_back_store,[],1);
  counts = hist(llr_back_large_store,centers);
  plot(centers,counts/sum(counts),'.-','color',colors(5,:));
  for l = 1:numel(prctiles),
    counts = hist(llr_fore_prctiles_store(l,:),centers);
    plot(centers,counts/sum(counts),'.-','color',colors(5+l,:));
  end
  legend('llr back','llr back edge','llr fore large','llr fore small','llr back large',prctilestrs{:});
  drawnow;
  
  
  fclose(fid);
  
end


%% some experiments we didn't train on

idxchosen = zeros(1,numel(experiments));
for i = 1:numel(experiments),
  [~,name] = fileparts(expdirs{i});
  idxchosen(i) = find(strcmp(name,expdirs_all));
end

experiments_notchosen = experiments_all(setdiff(1:numel(experiments_all),idxchosen));
ngal4_test = 4;
ncontrols_test = 4;
[experiments_test,ngal4_test_chosen,ncontrols_test_chosen] = choose_expdirs(experiments_notchosen,ngal4_test,ncontrols_test);

fid = fopen(expdirstestfilename,'w');
for i = 1:nexpdirs,
  [~,name] = fileparts(experiments_test(i).file_system_path);
  fprintf(fid,'%s\n',name);
end
fclose(fid);

fprintf('Run Ctrax on experiment %s.\n',experiments_test(1).file_system_path);
fprintf('Use old settings ann file, but change background model to read experiment background model from %s\n',outputFileName);
fprintf('Track long enough to save an ann file to %s\n',fullfile(ctrax_settingsdir,analysis_protocol,'settings.ann'));
fprintf('Run Ctrax on all experiments in file %s, look at results at some point.\n',expdirstestfilename);