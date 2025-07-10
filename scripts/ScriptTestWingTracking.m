% ScriptTestWingTracking

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlData';

%% list all experiments

data = SAGEListBowlExperiments('screen_type','primary','checkflags',true,'rootdir',rootdatadir);

%% select experiments from different lines, over different times

ngal4 = 100;
ncontrol = 10;

weight_order = {'genotype','date','rig','bowl'};
[experiments,~,~,idxchosen] = choose_expdirs(data,ngal4,ncontrol,'weight_order',weight_order);

%% output to file

filename = 'expdirs_testwingtracking_20120829.txt';
fid = fopen(filename,'w');
fprintf(fid,'%s\n',experiments.file_system_path);
fclose(fid);