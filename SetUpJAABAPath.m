% set up the JAABA path
% This shadows JAABA's version of this, b/c we don't want to add things already
% on the path, since that would unshadow some things (like ufmf_read_header())
% that we want to stay shadowed

fda_folder_path = fileparts(mfilename('fullpath')) ;
jaaba_folder_path = fullfile(fda_folder_path, 'JAABA') ;
jaaba_perframe_folder_path =  fullfile(jaaba_folder_path, 'perframe') ;
% Initialize all the paths.
addpath(jaaba_perframe_folder_path);  % in case we ever want to cd out of this dir
addpath(fullfile(jaaba_folder_path,'misc'));
%addpath(fullfile(jaaba_folder_path,'filehandling'));  % if we uncomment this, will
%unshadow ufmf_read_header()
addpath(fullfile(jaaba_perframe_folder_path,'larva_compute_perframe_features'));
addpath(fullfile(jaaba_perframe_folder_path,'compute_perframe_features'));
addpath(fullfile(jaaba_folder_path,'perframe','params'));
addpath(fullfile(jaaba_folder_path,'tests'));
st_dir = fullfile(jaaba_folder_path,'spaceTime'); 
addpath(genpath(st_dir));
adam_dir = fullfile(jaaba_folder_path,'spaceTime','adamMice'); 
rmpath(genpath(adam_dir));
