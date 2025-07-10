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
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';  
  rootdir_fixed = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/Fixed_AR/EP00005_rc8';
  %rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
end
analysis_protocol = '20110407';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};

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

%% experiments to analyze

experiment_params = struct;
experiment_params.experiment_name = 'FlyBowl_*';
% check for failures
experiment_params.checkflags = true;
% rearing condition
experiment_params.rearing_protocol = 'rc8';
% experimental protocol
experiment_params.experiment_protocol = 'EP0005.xml';
% root data dir
%experiment_params.rootdir = rootdir;
% manual_pf
experiment_params.manual_pf = {'P','U'};
% data types
experiment_params.data_type = {'registrationdata_circleCenterX',...
  'registrationdata_circleCenterY',...
  'registrationdata_circleRadius'};
experiment_params.dataset = 'score';
ngal4 = 25;
ncontrols = 25;
tmpparams = struct2paramscell(experiment_params);

% output file
expdirsfilename = fullfile(outputdir,expdirsFileStr);

%experiments_all = SAGEListBowlExperiments(tmpparams{:});
experiments_all = SAGEGetBowlData(tmpparams{:},'rootdir',rootdir);
[experiments_all.line] = deal(experiments_all.line_name);

%% sort experiments

ncontrols_all = nnz(strcmpi({experiments_all.line_name},'pBDPGAL4U'));
ngal4_all = numel(experiments_all) - ncontrols_all;
[experiments] = choose_expdirs(experiments_all,ngal4_all,ncontrols_all);

%% choose experiments with good background models, using ordering produced by sorting

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
    annfilename = fullfile(exp.file_system_path,dataloc_params.annfilestr);
    [bgcenter,bgstd,movie_height,movie_width] = ...
      read_ann(annfilename,'background_center','background_mad','movie_height','movie_width');
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
  
  if iscontrol,
    ncontrols_chosen = n;
  else
    ngal4_chosen = n;
  end
  
end

experiments_chosen = experiments(isgoodbg);

%% output chosen experiments to file

expdirs = {experiments_chosen.file_system_path};
nexpdirs = numel(expdirs);
fid = fopen(expdirsfilename,'w');
for i = 1:nexpdirs,
  fprintf(fid,'%s\n',expdirs{i});
end

%% command for learning

paramsFileStr = 'ExpBGFGModelParams.txt';
movieFileStr = dataloc_params.moviefilestr;
annFileStr = dataloc_params.annfilestr;
outputFileStr = 'ExpBGFGModelResults.pickle';
matFileStr = 'ExpBGFGModelResults.mat';

% file to write experiment directory names to
expdirsFileName = fullfile(outputdir,expdirsFileStr);

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
  annfilename = fullfile(expdir,dataloc_params.annfilestr);
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

[maxscore,i] = max(score);
LoG_thresh = edges(i+1);
fprintf('LoG_thresh = %f\n',LoG_thresh);

%% choose llr thresholds

expbgfgmodelfile = fullfile(settingsdir,analysis_protocol,dataloc_params.expbgfgmodelmatfilestr);
model = load(expbgfgmodelfile);

llr_back_store = [];
llr_fore_large_store = [];
llr_fore_small_store = [];
llr_back_edge_store = [];

lim = [-20,120];
nbins = 200;
nframessample = 20;
[edges,centers] = SelectHistEdges(nbins,lim,'linear');
[x,y] = meshgrid(1:movie_width,1:movie_height);
isnearedge = sqrt((x-movie_width/2).^2 + (y-movie_height/2).^2) >= .47*movie_width;

for expdiri = 1:nexpdirs,
  
  expdir = expdirs{expdiri};
  annfilename = fullfile(expdir,dataloc_params.annfilestr);
  moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
  trxfilename = fullfile(expdir,dataloc_params.ctraxfilestr);
  
  fprintf('%d: %s\n',expdiri,expdir);
  
  ann = read_ann(annfilename);
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
      llr_fore_large_store(end+1) = max(llrcurr);
      llr_fore_small_store(end+1) = min(llrcurr);
    end


    
  end
  
  hold off;
  counts = hist(llr_back_store,centers);
  plot(centers,counts/sum(counts),'k.-');
  hold on;
  counts = hist(llr_back_edge_store,centers);
  plot(centers,counts/sum(counts),'b.-');
  counts = hist(llr_fore_large_store,centers);
  plot(centers,counts/sum(counts),'r.-');
  counts = hist(llr_fore_small_store,centers);
  plot(centers,counts/sum(counts),'g.-');
  title(num2str(expdiri));
  drawnow;
  
  
  fclose(fid);
  
end

llr_back_large_store = max(llr_back_store,[],1);
counts = hist(llr_back_large_store,centers);
plot(centers,counts/sum(counts),'b.-');

%% some experiments we didn't train on

experiments_notchosen = experiments(~isgoodbg);
ngal4_test = 8;
ncontrols_test = 8;
for iscontrol = [false,true],
  idx = find(strcmpi({experiments_notchosen.line},'pBDPGAL4U') == iscontrol);
  if iscontrol,
    maxn = ncontrols_test;
  else
    maxn = ngal4_test;
  end
  for ii = 1:maxn,
    i = idx(ii);
    fprintf('%s\n',experiments_notchosen(i).file_system_path);
  end
end