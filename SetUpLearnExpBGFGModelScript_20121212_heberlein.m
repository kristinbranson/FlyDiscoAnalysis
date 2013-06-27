% set up for LearnExpBgFGModelScript.py

%% set up path, data locations

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
ctrax_settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax';
analysis_protocol = '20121212_non_olympiad_heberlein';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};
expdirsFileStr = 'expdirs.txt';
expfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20121212_non_olympiad_heberlein/expdirs_ChooseCtraxParameters_20121212.txt';

%% create directories

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

%% experiments to analyze

expdirs = importdata(expfile);
expdirs(cellfun(@isempty,expdirs)) = [];

% check that experiments exist
isbad = false(1,numel(expdirs));
for i = 1:numel(expdirs),
  if ~exist(expdirs{i},'dir'),
    isbad(i) = true;
    fprintf('Directory %s does not exist\n',expdirs{i});
  end
end

expdirs(isbad) = [];

%% read metadata

moviefilestr = dataloc_params.moviefilestr;
trxfilestr = dataloc_params.ctraxfilestr;
annfilestr = dataloc_params.annfilestr;

experiments_all = [];
issuccess = false(1,numel(expdirs));
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  res = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
  success = ~metadata.flag_aborted && exist(fullfile(expdir,moviefilestr),'file') && ...
    exist(fullfile(expdir,annfilestr),'file') && ...
    exist(fullfile(expdir,trxfilestr),'file');
  
  issuccess(i) = success;
  if success,
    res.file_system_path = expdir;
    if isempty(experiments_all),
      experiments_all = res;
    else
      experiments_all(end+1) = res; %#ok<SAGROW>
    end
  else
    fprintf('Excluding experiment %s from analysis, missing files or aborted\n',expdir);
  end
end

fprintf('Found %d / %d successful experiments\n',numel(experiments_all),numel(expdirs));

%% sort experiments

weight_order = {'rig','bowl','date','screen_type','screen_reason'};

[experiments,nchosen] = ChooseExperimentsCondition(experiments_all,numel(experiments_all),'weight_order',weight_order);

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);

%% choose experiments with good background models, using ordering produced by sorting

isgoodbg = false(1,numel(experiments));

maxn = 20;
n = 0;
for i = 1:nexpdirs,

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
  
experiments_chosen = experiments(isgoodbg);

%% output chosen experiments to file

expdirs = {experiments_chosen.file_system_path};
nexpdirs = numel(expdirs);
expdirsfilename = fullfile(outputdir,expdirsFileStr);
fid = fopen(expdirsfilename,'w');
for i = 1:nexpdirs,
  fprintf(fid,'%s\n',expdirs{i});
end
fclose(fid);

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

% LoG_thresh = 1.400000

%% choose llr thresholds

expbgfgmodelfile = fullfile(ctrax_settingsdir,analysis_protocol,'LearnCtraxParams',dataloc_params.expbgfgmodelmatfilestr);
model = load(expbgfgmodelfile);

% llr of pixels in the background image that are 
% not within always bkgd mask, above 90th percentile of llr
llr_back_store = [];
% largest llr within each fly
llr_fore_large_store = [];
% smallest llr within each fly
llr_fore_small_store = [];
% llr of pixels in the background image that are 
% not within always bkgd mask, near the image edge
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
      if isempty(llrcurr), continue; end
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
  xlabel('llr');
  ylabel('fraction');
  legend('bkgd','bkgd edge','fore large','fore small');
  drawnow;
  
  
  fclose(fid);
  
end

llr_back_large_store = max(llr_back_store,[],1);
counts = hist(llr_back_large_store,centers);
plot(centers,counts/sum(counts),'m.-');
legend('bkgd','bkgd edge','fore large','fore small','back large');

%% some experiments we didn't train on

experiments_notchosen = experiments(~isgoodbg);
maxn = 20;
for ii = 1:maxn,
  i = idx(ii);
  fprintf('%s\n',experiments_notchosen(i).file_system_path);
end
