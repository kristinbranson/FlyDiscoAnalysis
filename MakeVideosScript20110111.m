addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/misc;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/behavioralmicroarray;
addpath /groups/branson/bransonlab/projects/JCtrax/misc/;
addpath /groups/branson/bransonlab/projects/JCtrax/filehandling/;

%%

rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110111';
resultsdir = rootdir;
expdirs = {
  'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101220T142940'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101220T161020'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101220T175008'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlA_20101220T182340'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101220T142943'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101220T161026'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101220T175010'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlB_20101220T182342'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101220T142934'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101220T161017'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101220T175002'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlC_20101220T182334'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlD_20101220T161020'
  'pBDPGAL4U_TrpA_Rig1Plate01BowlD_20101220T182336'
}';


firstframe = 12000;
endframe = 13000;
nexpdirs = length(expdirs);
DEBUG = true;
moviefilestr = 'movie.ufmf';
trxfilestr = 'ctrax_results.mat';
annfilestr = 'movie.ufmf.ann';

%%

%avifiles = cell(1,nexpdirs);

for expdiri = expdiri:nexpdirs,

  expdir = fullfile(rootdir,expdirs{expdiri});

  % prepare to read movie
  moviefile = fullfile(expdir,moviefilestr);
  if ~exist(moviefile,'file'),
    error('Movie %s does not exist',moviefile);
  end

  % read annotation
  if DEBUG,
    tmp = regexp(expdir,'/?([^/]+)/?$','tokens','once');
    expdir_base = tmp{1};
    annfile = fullfile(resultsdir,expdir_base,annfilestr);
  else
    annfile = fullfile(expdir,annfilestr);
  end
  if ~exist(annfile,'file');
    warning('Annotation file %s does not exist',annfile);
    continue;
  end

  % % resize images read from annotation
  % annVarIsRead = false(1,size(var2ShowLabel,1));
  % fns = annVars;
  % for i = 1:length(fns),
  %   fn = fns{i};
  %   if isfield(ann,fn),
  %     ann.(fn) = reshape(ann.(fn),[nr,nc,numel(ann.(fn))/(nr*nc)]);
  %     annVarIsRead(i) = true;
  %   end
  % end
  % if ~isfield('background_center',ann),
  %   ann.background_center = ann.background_median;
  % end
  % if ~isfield('background_dev',ann),
  %   ann.background_dev = min(max(ann.background_mad,ann.bg_std_min),ann.bg_std_max);
  % end

  % read trx
  if DEBUG,
    trxfile = fullfile(resultsdir,expdir_base,trxfilestr);
  else
    trxfile = fullfile(expdir,trxfilestr);
  end
  if ~exist(trxfile,'file');
    warning('Mat file %s does not exist',trxfile);
    continue;
  end
  %trx = load_tracks(trxfile);
  
  avifiles{expdiri} = [expdir_base,sprintf('_frames%05dto%05d.avi',firstframe,endframe)];
  xvidfiles{expdiri} = fullfile(resultsdir,'figs',[expdir_base,sprintf('_frames%05dto%05d_xvid.avi',firstframe,endframe)]);

  [succeeded,avifiles{expdiri}] = ...
    make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avifiles{expdiri},'nzoomr',5,'nzoomc',3,...
    'boxradius',20,'taillength',25,'fps',30,'maxnframes',endframe-firstframe+1,'firstframe',firstframe);
  
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts bitrate=3000000',avifiles{expdiri},xvidfiles{expdiri});
  unix(cmd);
  
end
