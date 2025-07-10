
%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20111221/ctrax0.3_test';

%%


expdirs = importdata(fullfile(rootoutputdir,'expdirs.txt'));
fns_check = {'x','y','a','b','theta'};

for i = 1:numel(expdirs),
  
  fprintf('%s...\n',expdirs{i});
  trx = load_tracks(fullfile(rootoutputdir,expdirs{i},'test_ctrax_results.mat'));
  for j = 1:numel(fns_check),
    if any(isinf([trx.(fns_check{j})])),
      fprintf('%s %s has infs\n',expdirs{i},fns_check{j});
    end
    if any(isnan([trx.(fns_check{j})])),
      fprintf('%s %s has nans\n',expdirs{i},fns_check{j});
    end
    
  end
  
end
