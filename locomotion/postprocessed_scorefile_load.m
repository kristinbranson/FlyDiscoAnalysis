function [scores] = postprocessed_scorefile_load(expdir,scorefilename)

load(fullfile(expdir,scorefilename));

for fly = 1:numel(allScores.postprocessed)
scores{fly} = allScores.postprocessed{fly}(allScores.tStart(fly):allScores.tEnd(fly));
end