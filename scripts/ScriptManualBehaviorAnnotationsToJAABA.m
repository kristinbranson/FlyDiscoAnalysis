
rootdir = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/ManualFlyTrackerToJAABA';
% for per-frame parameters
perframeparamsfile = '/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/settings/current_non_olympiad_dickson_VNC/perframe_params.txt';
analysis_protocol = 'current_non_olympiad_dickson_VNC';
basejabfile = fullfile(rootdir,'manual.jab');
baseexpdirs = mydir(rootdir,'isdir',true);
expnames = cell(size(baseexpdirs));
expdirs = cell(size(baseexpdirs));
nexps = numel(expdirs);
movieext = '.mp4';
moviestr = ['movie',movieext];
manualfiles = cell(size(expnames));
trxfilestr = 'trx.mat';
perframedirstr = 'perframe';
outtrxfilestr = 'converted_trx.mat'; 
load('/groups/branson/home/bransonk/behavioranalysis/code/FlyDiscoAnalysis/ManualFlyTrackerToJAABA/JAABAFeatures.mat','windowFeaturesParams');


for i = 1:numel(baseexpdirs),
  [~,expnames{i}] = fileparts(baseexpdirs{i});
  expdirs{i} = fullfile(baseexpdirs{i},[expnames{i},'_JAABA']);
  assert(exist(expdirs{i},'dir')>0);
  moviefile0 = fullfile(baseexpdirs{i},[expnames{i},movieext]);
  assert(exist(moviefile0,'file')>0);
  moviefile = fullfile(expdirs{i},moviestr);
  if ~exist(moviefile,'file'),
    unix(sprintf('ln -s %s %s',moviefile0,moviefile));
  end
  assert(exist(moviefile,'file')>0);
  manualfiles{i} = fullfile(baseexpdirs{i},[expnames{i},'-actions.mat']);
end

perframe_params = ReadParams(perframeparamsfile);

% fix issues with missing files
forcecompute = false;

fns = {
  'timestamps'
  'dt'
  'x'
  'y'
  'theta'
  'a'
  'b'
  'xwingl'
  'ywingl'
  'xwingr'
  'ywingr'
  'x_mm'
  'y_mm'
  'a_mm'
  'b_mm'
  'theta_mm'
  };

for expi = 1:nexps,
  outtrxfile = fullfile(expdirs{expi},outtrxfilestr);
  perframedir = fullfile(expdirs{expi},perframedirstr);

  if exist(outtrxfile,'file'),
    if forcecompute,
      delete(outtrxfile); 
    else
      continue; 
    end
  end
  
  ftfile = fullfile(baseexpdirs{expi},[expnames{expi},'-track.mat']);
  trxfile = fullfile(expdirs{expi},trxfilestr);
  % row arrays
  td = load(trxfile);
  for i = 1:numel(td.trx),
    for j = 1:numel(fns),
      td.trx(i).(fns{j}) = td.trx(i).(fns{j})(:)';
    end
  end
  
  if ~exist(perframedir,'dir'),
    mkdir(perframedir);
  end
  calibfile = fullfile(baseexpdirs{expi},[expnames{expi},'-calibration.mat']);
  cald = load(calibfile);
  
  arena = struct;
  arena.x = cald.calib.centroids(2);
  arena.y = cald.calib.centroids(1);
  arena.r = cald.calib.r;
  
  if exist(outtrxfile,'file'),
    delete(outtrxfile);
  end

  outtrx = FlyTracker2WingTracking_helper(ftfile,td,perframedir,outtrxfile,perframe_params,arena,[]);
end


expi = 1;
md = load(manualfiles{expi});
behs = md.behs;
behs = strrep(behs,'-','_');
nbehs = numel(behs);
jabfiles = cell(1,nbehs);

for behi = 1:nbehs,
  jabfiles{behi} = fullfile(rootdir,[behs{behi},'.jab']);
end
for behi = 1:nbehs,
  if exist(jabfiles{behi},'file'),
    if forcecompute,
      delete(jabfiles{behi});
    else
      continue;
    end
  end
  jd = JLabelData();
  jd.isInteractive = false;
  jd.openJabFile(basejabfile,false); % false means not groundtruthing mode
  jd.renameBehavior(jd.labelnames{1},behs{behi}); % change behavior name
  
  for expi = 1:nexps,
    jd.AddExpDir(expdirs{expi});
  end
  
  jd.saveJabFile(jabfiles{behi});
  
end

trainedjabfiles = cell(1,nbehs);
for behi = 1:nbehs,
  [p,n,ext] = fileparts(jabfiles{behi});
  trainedjabfiles{behi} = fullfile(p,[n,'_trained',ext]);
end

maxnsamples = 10000;
nflies = 2;

label_pred_counts = nan([nflies,nexps,2,2,nbehs]);
fpr = nan([1,nbehs]);
fnr = nan([1,nbehs]);
fpr_per_exp = nan(nexps,nbehs);
fnr_per_exp = nan(nexps,nbehs);
fpr_per_fly = nan(nflies,nexps,nbehs);
fnr_per_fly = nan(nflies,nexps,nbehs);
nframespos_all = nan([nflies,nexps,nbehs]);
nframesneg_all = nan([nflies,nexps,nbehs]);

labelvals = [2,1];

for behi = 2:nbehs,
  
  jd = JLabelData();
  jd.isInteractive = false;
  jd.openJabFile(jabfiles{behi},false); % false means not groundtruthing mode
  
  nframespos = zeros(nflies,nexps);
  nframesneg = zeros(nflies,nexps);
  for expi = 1:nexps,
    md = load(manualfiles{expi});
    nframes = jd.endframes_per_exp{expi}(1);
    for flies = 1:jd.nflies_per_exp(expi),
      labelidx = FlyTrackerBouts2LabelIdx(md.bouts(flies,behi),nframes);
      nframespos(flies,expi) = nnz(labelidx==1);
      nframesneg(flies,expi) = nnz(labelidx==2);
    end
  end
  nframespos_all(:,:,behi) = nframespos;
  nframesneg_all(:,:,behi) = nframesneg;
  if sum(nframespos(:)) == 0,
    warning('No positive labels, skipping %s\n',behs{behi});
    continue;
  end
  if sum(nframesneg(:)) == 0,
    warning('No negative labels, skipping %s\n',behs{behi});
    continue;
  end  
  
  maxnsamplesperexp = maxnsamples / nexps;
  maxnsamplesperfly = maxnsamplesperexp / nflies;
  
  if sum(nframespos(:)) <= maxnsamples,
    psample_pos = ones(nflies,nexps);
  else
    wfly = 1./(max(1,nframespos));
    wfly(nframespos==0) = 0;
    wfly = wfly ./ sum(wfly(:)) * numel(nframespos);
    psamples = permute((1:maxnsamples)/maxnsamples,[1,3,2]);
    expnsample_c = sum(sum(min(1,psamples.*wfly).*nframespos,1),2);
    assert(min(expnsample_c)<maxnsamples && max(expnsample_c)>maxnsamples);
    ci = find(expnsample_c<=maxnsamples,1,'last');
    psample_pos = min(1,psamples(ci).*wfly);
    %expnsample = sum(sum(psample_pos.*nframespos));
  end
  
  if sum(nframesneg(:)) <= maxnsamples,
    psample_neg = ones(nflies,nexps);
  else
    wfly = 1./(max(1,nframesneg));
    wfly(nframesneg==0) = 0;
    wfly = wfly ./ sum(wfly(:)) * numel(nframesneg);
    expnsample_c = sum(sum(min(1,psamples.*wfly).*nframesneg,1),2);
    ci = find(expnsample_c<=maxnsamples,1,'last');
    assert(~isempty(ci));
    psample_neg = min(1,psamples(ci).*wfly);
  end
  
  for expi = 1:nexps,
    md = load(manualfiles{expi});
    
    assert(strcmp(jd.expdirs{expi},expdirs{expi}));
    assert(all(jd.firstframes_per_exp{expi}==1));
    nframes = jd.endframes_per_exp{expi}(1);
    assert(all(jd.endframes_per_exp{expi}==nframes));
    ts = 1:nframes;
    
    for flies = 1:jd.nflies_per_exp(expi),
      labelidx = FlyTrackerBouts2LabelIdx(md.bouts(flies,behi),nframes);
      dosample = rand(1,nframes);
      dosample(labelidx==1) = dosample(labelidx==1) <= psample_pos(flies,expi);
      dosample(labelidx==2) = dosample(labelidx==2) <= psample_neg(flies,expi);
      labelidx(~dosample) = 0;
      
      jd.setLabelsRaw(expi,flies,ts,1,labelidx,1);
    end
  end
  
  jd.setCurrentTarget(1,1);
  jd.setWindowFeaturesParams(windowFeaturesParams);
  jd.saveJabFile(trainedjabfiles{behi});
  jd.Train();
  jd.saveJabFile(trainedjabfiles{behi});
  
  [classifierinfo,allScores] = JAABADetect(expdirs,'jabfiles',{trainedjabfiles{behi}});
  counts = zeros(nflies,nexps,2,2); % fly, exp, label, pred
  for expi = 1:nexps,
    md = load(manualfiles{expi});
    for flies = 1:jd.nflies_per_exp(expi),
      labelidx = FlyTrackerBouts2LabelIdx(md.bouts(flies,behi),nframes);
      for labeli = 1:2,
        label = labelvals(labeli);
        for pred = 0:1,
          counts(flies,expi,labeli,pred+1) = nnz(labelidx==label & allScores{expi}{1}.postprocessed{flies}==pred);
        end
      end
    end
  end
  %                 FP: label=0,pred=1     N: label=0
  fpr_per_fly_curr = counts(:,:,1,2)./sum(counts(:,:,1,:),4);
  fpr_per_exp_curr = sum(counts(:,:,1,2),1)./sum(sum(counts(:,:,1,:),4),1);
  fpr_curr = sum(sum(counts(:,:,1,2),1),2)./sum(sum(sum(counts(:,:,1,:),4),1),2);
  %                  FN: label=1,pred=0     N: label=1
  fnr_per_fly_curr = counts(:,:,2,1)./sum(counts(:,:,2,:),4);
  fnr_per_exp_curr = sum(counts(:,:,2,1),1)./sum(sum(counts(:,:,2,:),4),1);
  fnr_curr = sum(sum(counts(:,:,2,1),1),2)./sum(sum(sum(counts(:,:,2,:),4),1),2);
  
  label_pred_counts(:,:,:,:,behi) = counts;
  fpr(behi) = fpr_curr;
  fnr(behi) = fnr_curr;
  fpr_per_exp(:,behi) = fpr_per_exp_curr;
  fnr_per_exp(:,behi) = fnr_per_exp_curr;
  fpr_per_fly(:,:,behi) = fpr_per_fly_curr;
  fnr_per_fly(:,:,behi) = fnr_per_fly_curr;
  
  fprintf('%s: FPR = %f, FNR = %f\n',behs{behi},fpr_curr,fnr_curr);
  fprintf('FPR per fly: %s\n',mat2str(fpr_per_fly_curr,2));
  fprintf('FNR per fly: %s\n',mat2str(fnr_per_fly_curr,2));
  
end

for behi = 1:nbehs,
  fprintf('\n%s:\nFPR = %f, FNR = %f\n',behs{behi},fpr(behi),fnr(behi));
  fprintf('N pos labels = %d\n',sum(sum(nframespos_all(:,:,behi))));
  fprintf('N neg labels = %d\n',sum(sum(nframesneg_all(:,:,behi))));
  fprintf('N pos labels per fly = %s\n',mat2str(nframespos_all(:,:,behi)));
  fprintf('N neg labels per fly = %s\n',mat2str(nframesneg_all(:,:,behi)));
  fprintf('FPR per fly: %s\n',mat2str(fpr_per_fly(:,:,behi),2));
  fprintf('FNR per fly: %s\n',mat2str(fnr_per_fly(:,:,behi),2));
  fprintf('FPR per exp: %s\n',mat2str(fpr_per_exp(:,behi),2));
  fprintf('FNR per exp: %s\n',mat2str(fnr_per_exp(:,behi),2));
end


