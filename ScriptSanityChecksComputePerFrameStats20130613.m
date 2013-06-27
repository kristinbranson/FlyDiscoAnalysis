%% check not any behavior
scoresfiles = {
  'scoresBackup.mat'
  'scoresTouch.mat'
  'scoresWingGrooming.mat'
  'scores_AttemptedCopulation.mat'
  'scores_BodyTurns.mat'
  'scores_Chasev7.mat'
  'scores_Copulation.mat'
  'scores_Crabwalk2.mat'
  'scores_Crabwalk3.mat'
  'scores_Jump.mat'
  'scores_Righting.mat'
  'scores_Stops.mat'
  'scores_Walk.mat'
  'scores_WingExtension.mat'
  'scores_pivot_center.mat'
  'scores_pivot_tail.mat'
  'scores_wingflick.mat'};
nbehaviors = numel(scoresfiles);

td = load(fullfile(expdir,'registered_trx.mat'));
nflies = numel(td.trx);
behaviorsdetected = cell(1,nflies);
for fly = 1:nflies,
  behaviorsdetected{fly} = false(nbehaviors,td.trx(fly).nframes);
end
for i = 1:nbehaviors,
  sc = load(fullfile(expdir,scoresfiles{i}));
  for fly = 1:nflies,
    tstart = td.trx(fly).firstframe;
    tend = min(numel(sc.allScores.postprocessed{fly}),td.trx(fly).endframe);
    behaviorsdetected{fly}(i,1:tend-tstart+1) = sc.allScores.postprocessed{fly}(tstart:tend);
  end
end

% % make sure that postprocessed is the same as t0s and t1s
% for i = 1:nbehaviors,
%   sc = load(fullfile(expdir,scoresfiles{i}));
%   for fly = 1:nflies,
%     tstart = td.trx(fly).firstframe;
%     tend = min(numel(sc.allScores.postprocessed{fly}),td.trx(fly).endframe);
%     tmp1 = sc.allScores.postprocessed{fly}(tstart:tend);
%     tmp2 = false(1,td.trx(fly).nframes);
%     for j = 1:numel(sc.allScores.t0s{fly}),
%       i0 = sc.allScores.t0s{fly}(j)-td.trx(fly).firstframe+1;
%       i1 = sc.allScores.t1s{fly}(j)-td.trx(fly).firstframe+1;
%       tmp2(i0:i1-1) = true;
%     end
%     if any(tmp1~=tmp2),
%       fprintf('Difference for behavior %s, fly %d\n',scoresfiles{i},fly);
%     end
%   end
% end

nonedetected = cell(1,nflies);
for fly = 1:nflies,
  nonedetected{fly} = ~any(behaviorsdetected{fly},1);
end

fractimenotanybehavior_perfly = cellfun(@nnz,nonedetected)./cellfun(@numel,nonedetected);
fractimenotanybehavior_perexp = mean(fractimenotanybehavior_perfly);

sd = load(fullfile(expdir,'stats_perframe.mat'));
fprintf('Comparison to not any behavior (per fly):\n');
fractimenotanybehavior_perfly - sd.statsperfly.fractime_flyany_framenotanybehavior.mean

fprintf('Comparison to not any behavior (per exp):\n');
fractimenotanybehavior_perexp - sd.statsperexp.fractime_flyany_framenotanybehavior.meanmean

%% check chase not wing extension

chasei = find(strcmp(scoresfiles,'scores_Chasev7.mat'));
wingextensioni = find(strcmp(scoresfiles,'scores_WingExtension.mat'));

chasenotwe = cell(1,nflies);
chasenotwe_male = cell(1,nflies);
scc = load(fullfile(expdir,scoresfiles{chasei}));
scw = load(fullfile(expdir,scoresfiles{wingextensioni}));
nmale = zeros(1,nflies);
for fly = 1:nflies,
  ismale = strcmp(td.trx(fly).sex,'M');
  tstart = td.trx(fly).firstframe;
  tend = min(numel(sc.allScores.postprocessed{fly}),td.trx(fly).endframe);
  chasenotwe{fly} = scc.allScores.postprocessed{fly}(tstart:tend) & ~scw.allScores.postprocessed{fly}(tstart:tend);
  chasenotwe_male{fly} = chasenotwe{fly}(ismale);
  nmale(fly) = nnz(ismale);
end
fractimechasenotwe_perfly = cellfun(@nnz,chasenotwe)./cellfun(@numel,chasenotwe);
fractimechasenotwe_perexp = mean(fractimechasenotwe_perfly);

fractimechasenotwe_male_perfly = cellfun(@nnz,chasenotwe_male)./cellfun(@numel,chasenotwe_male);
fractimechasenotwe_male_perexp = weighted_mean(fractimechasenotwe_male_perfly',nmale');


fprintf('Comparison to chase not wing extension (per fly):\n');
fractimechasenotwe_perfly - sd.statsperfly.fractime_flyany_framechase_notwingextension.mean

fprintf('Comparison to chase not wing extension (per exp):\n');
fractimechasenotwe_perexp - sd.statsperexp.fractime_flyany_framechase_notwingextension.meanmean

fprintf('Comparison to male chase not wing extension (per fly):\n');
tmp = fractimechasenotwe_male_perfly - sd.statsperfly.fractime_flymale_framechase_notwingextension.mean;
tmp(isnan(fractimechasenotwe_male_perfly) & isnan(sd.statsperfly.fractime_flymale_framechase_notwingextension.mean)) = 0

fprintf('Comparison to male chase not wing extension (per exp):\n');
tmp = fractimechasenotwe_male_perexp - sd.statsperexp.fractime_flymale_framechase_notwingextension.meanmean;
tmp(isnan(fractimechasenotwe_male_perexp) & isnan(sd.statsperexp.fractime_flymale_framechase_notwingextension.meanmean)) = 0

%% check duration of backups

backi = find(strcmp(scoresfiles,'scoresBackup.mat'));
sc = load(fullfile(expdir,scoresfiles{backi}));
durs = cell(1,nflies);
for fly = 1:nflies,
  tstart = td.trx(fly).firstframe;
  tend = min(numel(sc.allScores.postprocessed{fly}),td.trx(fly).endframe);
  [t0s,t1s] = get_interval_ends(sc.allScores.postprocessed{fly}(tstart:tend));
  durs{fly} = td.trx(fly).timestamps(t1s)-td.trx(fly).timestamps(t0s);
end

meandur_perfly = cellfun(@mean,durs);
fprintf('Comparison between dur backup:\n');
meandur_perfly - sd.statsperfly.duration_flyany_framebackup.mean

%%


% 
% {  'labelsBackup'
%   'labelsTouch'
%   'labelsWingGrooming'
%   'labels_AttemptedCopulation'
%   'labels_BodyTurns'
%   'labels_Chasev7'
%   'labels_Copulation'
%   'labels_Crabwalk2'
%   'labels_Crabwalk3'
%   'labels_Jump'
%   'labels_Righting'
%   'labels_Stops'
%   'labels_Walk'
%   'labels_WingExtension'
%   'labels_pivot_center'
%   'labels_pivot_tail'
%   'labels_wingflick'};
