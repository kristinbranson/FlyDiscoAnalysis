%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% data locations

if ispc,
  
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  analysis_protocol = '20110202_pc';
  dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,'dataloc_params.txt'));
  expdir = 'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734';
  expdir_read = fullfile(dataloc_params.rootreaddir,expdir);
  
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  protocol = '';
  analysis_protocol = '20110211';
  expdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
end
% [readframe,nframes,fid,vidinfo] = get_readframe_fcn(fullfile(expdir_read,'movie.ufmf'));
% fclose(fid);

%% run

FlyBowlMakeCtraxResultsMovie(expdir);