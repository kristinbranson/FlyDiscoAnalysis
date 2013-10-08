%% set up path

% note that this script only fixes fractime_flyany_framenotanybehavior, not
% other framenotanybehavior statistics

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130920.mat';
%outdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130920.mat';

%% parameters

scoresfilestrs = {
%   'scores_AttemptedCopulation.mat'
%   'scoresBackup.mat'
%   'scores_Chasev7.mat'
%   'scores_Copulation.mat'
%   'scores_Crabwalk2.mat'
  'scores_Jump.mat'
%   'scores_pivot_center.mat'
%   'scores_pivot_tail.mat'
%   'scores_Righting.mat'
%   'scores_Walk.mat'
%   'scores_Stops.mat'
%   'scoresTouch.mat'
  'scores_WingExtension.mat'
%   'scoresWingGrooming.mat'
%   'scores_wingflick.mat'
  };
% 
% scoresfilestrs = {
% 'scores_AttemptedCopulation.mat'
% 'scoresBackup.mat'
% 'scores_BodyTurns.mat'
% 'scores_Chasev7.mat'
% 'scores_Copulation.mat'
% 'scores_Crabwalk3.mat'
% 'scores_Crabwalk2.mat'
% 'scores_Jump.mat'
% 'scores_pivot_center.mat'
% 'scores_pivot_tail.mat'
% 'scores_Righting.mat'
% 'scores_Walk.mat'
% 'scores_Stops.mat'
% 'scoresTouch.mat'
% 'scores_WingExtension.mat'
% 'scoresWingGrooming.mat'
% 'scores_wingflick.mat'
% }

minnframes = 500;
maxconfidence = 0.25;
nbehaviors = numel(scoresfilestrs);
minnexpsperset = 2;

%% load in data

load(datafile,'metadata','linestats');

%% compute 

nexps = numel(metadata);
isproblem = false(nbehaviors,nexps);

nlowconfidencepos = zeros(nbehaviors,nexps);
nlowconfidenceneg = zeros(nbehaviors,nexps);
npos = zeros(nbehaviors,nexps);
nneg = zeros(nbehaviors,nexps);
for expi = 1:nexps,
  
  expdir = metadata(expi).file_system_path;
  for j = 1:numel(scoresfilestrs),
    filename = fullfile(expdir,scoresfilestrs{j});
    if ~exist(filename,'file'),
      warning('File %s does not exist.',filename);
      isproblem(j,expi) = true;
      continue;
    end
    sd = load(filename);
    nflies = numel(sd.allScores.tStart);
    nframes = sd.allScores.tEnd-sd.allScores.tStart+1;
    for fly = 1:nflies,
      if nframes(fly) < minnframes,
        continue;
      end
      npos(j,expi) = npos(j,expi) + nnz(sd.allScores.postprocessed{fly}==1);
      nneg(j,expi) = nneg(j,expi) + nnz(sd.allScores.postprocessed{fly}==0);
      nlowconfidencepos(j,expi) = nlowconfidencepos(j,expi) + ...
        nnz(sd.allScores.scores{fly}(sd.allScores.postprocessed{fly}==1)/sd.allScores.scoreNorm < maxconfidence);
      nlowconfidenceneg(j,expi) = nlowconfidenceneg(j,expi) + ...
        nnz(sd.allScores.scores{fly}(sd.allScores.postprocessed{fly}==0)/sd.allScores.scoreNorm > -maxconfidence);
    end
  
    fprintf('Experiment %d / %d, %s, fracpos = %f, fracneg = %f\n',...
      expi,nexps,scoresfilestrs{j},nlowconfidencepos(j,expi)/npos(j,expi),...
      nlowconfidenceneg(j,expi)/nneg(j,expi));

  end
  
end

lowconfidencerate = (nlowconfidencepos+nlowconfidenceneg)./(npos+nneg);
posrate = npos./(npos+nneg);

%% combine to get per-set values

[set_names,firstidx,setidx] = unique({metadata.set});
nsets = numel(set_names);
posrate_set = nan(nbehaviors,nsets);
lowconfidencerate_set = nan(nbehaviors,nsets);
isokset = false(1,nsets);
for seti = 1:nsets,
  posrate_set(:,seti) = mean(posrate(:,setidx==seti),2);
  lowconfidencerate_set(:,seti) = mean(lowconfidencerate(:,setidx==seti),2);
  isokset(seti) = nnz(setidx==seti) >= minnexpsperset;
end

%% combine to get per-line values

setmetadata = metadata(firstidx);
[~,set2lineidx] = ismember({setmetadata.line_name},linestats.line_names);
nlines = numel(linestats.line_names);
posrate_line = nan(nbehaviors,nlines);
lowconfidencerate_line = nan(nbehaviors,nlines);
for linei = 1:numel(linestats.line_names),
  idx = set2lineidx==linei & isokset;
  posrate_line(:,linei) = mean(posrate_set(:,idx),2);
  lowconfidencerate_line(:,linei) = mean(lowconfidencerate_set(:,idx),2);
end

medianlowconfidencerate = nanmedian(lowconfidencerate_line,2);
normlowconfidencerate_line = zeros(nbehaviors,nlines);
orderedlines = nan(nbehaviors,nlines);
for i = 1:nbehaviors,
  idx = lowconfidencerate_line(i,:) >= medianlowconfidencerate(i);
  normlowconfidencerate_line(i,idx) = lowconfidencerate_line(i,idx)./posrate_line(i,idx);
  [~,orderedlines(i,:)] = sort(normlowconfidencerate_line(i,:),'descend');
end

for i = 1:nbehaviors,
  
  fprintf('\n%s:\n',scoresfilestrs{i});
  for jj = 1:2,
    j = orderedlines(i,jj);
    idx = find(strcmp({metadata.line_name},linestats.line_names{j}));
    tmp1 = mean(lowconfidencerate(i,idx));
    tmp2 = mean(posrate(i,idx));
    fprintf('%s: low confidence rate = %f (should be ~ %f), pos rate = %f (should be ~ %f)\n',...
      linestats.line_names{j},lowconfidencerate_line(i,j),tmp1,posrate_line(i,j),tmp2);    
  end  
  
end

% scores_Jump.mat:
% GMR_23A07_AE_01: low confidence rate = 0.008204 (should be ~ 0.008146), pos rate = 0.000363 (should be ~ 0.000358)
% GMR_54A12_AE_01: low confidence rate = 0.002391 (should be ~ 0.002391), pos rate = 0.000129 (should be ~ 0.000129)
% 
% scores_WingExtension.mat:
% GMR_57G11_AD_01: low confidence rate = 0.001395 (should be ~ 0.001050), pos rate = 0.000000 (should be ~ 0.000000)
% GMR_55F04_AE_01: low confidence rate = 0.001542 (should be ~ 0.001426), pos rate = 0.000007 (should be ~ 0.000008)


%% save results

save -v7.3 FracTimeLowConfidence20130924.mat lowconfidencerate* posrate* medianlowconfidencerate normlowconfidencerate_line orderedlines metadata isokset
