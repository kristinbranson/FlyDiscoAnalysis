
%% set up paths

addpath(genpath('/groups/branson/bransonlab/projects/olympiad/jfrc_metadata_tools/jdbc'));
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
%sage_params_path = '/groups/branson/home/bransonk/SAGEConfig/SAGEReadParams.txt';
sage_params_path = 'prod.params';

%% get data

infilename = 'PropagateMetadataExperimentNames_Small.txt';
outfilename = 'PropagateMetadataTODO.txt';
params = {'debug',0,'sage_params_path',sage_params_path};

propagateMetadataFixes(infilename,outfilename,params{:});

