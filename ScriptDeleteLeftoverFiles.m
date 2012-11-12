% delete files from screen data that shouldn't be there

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;


filelist = 'filesownedbybransonk.txt';

%% find all base names

files = importdata(filelist);
basenames = cell(1,numel(files));
for i = 1:numel(files),
  [~,basenames{i}] = myfileparts(files{i});
end
[unique_basenames,~,nameidx] = unique(basenames);

fprintf('%s\n',unique_basenames{:});

%% choose which ones to delete

names_delete = {
  'Thumbs.db'
  'bak_110804-145642_testmovie.ufmf.ann'
  'bak_110804-150033_testmovie.ufmf.ann'
  'bak_110804-152849_testmovie.ufmf.ann'
  'bak_110804-152858_testmovie.ufmf.ann'
  'bak_110804-153134_testmovie.ufmf.ann'
  'bak_110804-153311_testmovie.ufmf.ann'
  'bak_110804-153846_testmovie.ufmf.ann'
  'bak_110804-154404_testmovie.ufmf.ann'
  'bak_110804-154532_testmovie.ufmf.ann'
  'bak_110804-154724_testmovie.ufmf.ann'
  'bak_110804-154836_testmovie.ufmf.ann'
  'bak_110804-154945_testmovie.ufmf.ann'
  'bak_110804-155100_testmovie.ufmf.ann'
  'bak_110804-203927_testmovie.ufmf.ann'
  'bak_110804-204121_testmovie.ufmf.ann'
  'bak_110804-204152_testmovie.ufmf.ann'
  'bak_110804-204227_testmovie.ufmf.ann'
  'bak_110804-204405_testmovie.ufmf.ann'
  'bak_110804-204416_testmovie.ufmf.ann'
  'bak_110804-204819_testmovie.ufmf.ann'
  'bak_110804-204835_testmovie.ufmf.ann'
  'bak_110805-095100_testmovie.ufmf.ann'
  'bak_110808-092639_testmovie.ufmf.ann'
  'ctrax_results_movie_pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110923T142257.avi.tmp'
  'ctrax_results_movie_pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110923T142257_small.avi'
  'ctrax_results_movie_pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110923T142257_temp.avi'
  'movie_ctraxdiagnostics.txt'
  'new_ctrax_diagnostics.txt'
  'new_ctrax_results.mat'
  'newmovie.ufmf.ann'
  'stats_perframe.txt.new'
  'subtitles.srt'
  'test_ctrax_diagnostics.txt'
  'test_ctrax_results.mat'
  'testmovie.ufmf.ann'
  };

dodelete = ismember(basenames,names_delete);
fprintf('Deleting %d files, keeping %d files\n',nnz(dodelete),nnz(~dodelete));

isproblem = false(1,numel(files));
for i = 1:numel(files),
  cmd = sprintf('chmod g+w %s',files{i});
  [status,result] = unix(cmd);
  if status,
    fprintf('Error changing write permissions for file %s: %s\n',files{i},result);
    isproblem(i) = true;
  end
end
% for i = 1:numel(files),
%   cmd = sprintf('chmod o-w %s',files{i});
%   [status,result] = unix(cmd);
%   if status,
%     fprintf('Error changing write permissions for file %s: %s\n',files{i},result);
%   end
% end

fid = fopen('filesownedbybransonk_delete.txt','w');
fprintf(fid,'%s\n',files{dodelete});
fclose(fid);

fid = fopen('filesownedbybransonk_chown.txt','w');
fprintf(fid,'%s\n',files{~dodelete});
fclose(fid);
