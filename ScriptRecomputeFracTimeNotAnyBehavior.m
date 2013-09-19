%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130912.mat';

%% parameters

scoresfilestrs = {
  'scores_AttemptedCopulation.mat'
  'scoresBackup.mat'
  'scores_Chasev7.mat'
  'scores_Copulation.mat'
  'scores_Crabwalk2.mat'
  'scores_Jump.mat'
  'scores_pivot_center.mat'
  'scores_pivot_tail.mat'
  'scores_Righting.mat'
  'scores_Walk.mat'
  'scores_Stops.mat'
  'scoresTouch.mat'
  'scores_WingExtension.mat'
  'scoresWingGrooming.mat'
  'scores_wingflick.mat'
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

%% load in data

load(datafile);

%% recompute fractime_flyany_framenotanybehavior

nexps = numel(metadata);
isproblem = false(1,nexps);
fractimenotanybehavior = nan(1,nexps);
for expi = 1:nexps,
  
  fprintf('Experiment %d / %d\n',expi,nexps);
  
  expdir = metadata(expi).file_system_path;
  for j = 1:numel(scoresfilestrs),
    filename = fullfile(expdir,scoresfilestrs{j});
    if ~exist(filename,'file'),
      warning('File %s does not exist.',filename);
      isproblem(expi) = true;
      break;
    end
    sd = load(filename);
    if j == 1,
      tStarts = sd.allScores.tStart;
      tEnds = sd.allScores.tEnd;
      nframes = tEnds-tStarts+1;
      nflies = numel(tStarts);
      isbehavior = cell(1,nflies);
      for fly = 1:nflies,
        isbehavior{fly} = false(1,nframes(fly));
      end
    else
      if numel(tStarts) ~= numel(sd.allScores.tStart) || ...
          numel(tEnds) ~= numel(sd.allScores.tEnd),
        warning('Scores file %s does not match previous scores files.\n',filename);
        isproblem(expi) = true;
        break;
      end
    end

    for fly = 1:nflies,
      if nframes(fly) >= minnframes,
        isbehavior{fly} = isbehavior{fly} | sd.allScores.postprocessed{fly}(tStarts(fly):tEnds(fly))>0;
      end
    end
  end
  
  if isproblem(expi),
    continue;
  end
  
  isbehavior = [isbehavior{:}];
  fractimenotanybehavior(expi) = nnz(~isbehavior) / numel(isbehavior);
  
end