%% set up paths

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

%% data locations

if ispc,
else
  protocol = 'CtraxTest20110202';
  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol',protocol,'subreadfiles',{'ctrax_results.mat'});
  expdir = expdirs{1};
  [readframe,nframes,fid,vidinfo] = get_readframe_fcn(fullfile(expdir_reads{1},'movie.ufmf'));
  fclose(fid);
end

%%
obj = Trx();
obj.AddExpDir(expdir,vidinfo);