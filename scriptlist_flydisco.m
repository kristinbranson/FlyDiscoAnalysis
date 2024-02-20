% scripts for FlyDisco Data
% AR 20210316
modpath
%% review flytracker output for NaN gaps
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654/flytracker/movie/movie-track.mat';
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210408_testingpipeline/VNC_YNA_K_162984_RigD_20210405T152704/movie-track.mat';
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/nochr_TrpA71G01_Unknown_RigA_20201216T162938/movie-trk.mat';
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/nochr_TrpA71G01_Unknown_RigA_20201216T162938/movie-track.mat'
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/nochr_TrpA65F12_Unknown_RigD_20201216T175902/movie-track.mat'
% filename = '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS25523_RigB_20210324T144555/movie-track.mat'
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/VNC_JRC_SS29892_RigD_20210323T134743/movie-track.mat';
% filename = 'VNC_JRC_SS32371_RigC_20210325T141559'
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/VNC_JHS_K_85321_RigB_20210325T152834/movie-track.mat'; %none
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/VNC_JRC_SS29662_RigA_20210324T133352/movie-track.mat'; % none
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/VNC_EXT_VGLUT-GAL4_RigC_20210324T161114/movie-track.mat'; % none
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/VNC_JRC_SS33430_RigB_20210325T133338/movie-track.mat'; %none
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/cx_GMR_OL0046B_CsChr_RigA_20151020T135455/movie-track.mat'; %none
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654/movie-track.mat'; %already used - track a different social movie, 
% filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210518_testingFlyTrackerMods/nochr_TrpA71G01_Unknown_RigB_20201216T163123/movie-track.mat';
filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210518_testingFlyTrackerMods/VNC_JHS_K_85321_RigA_20210322T164447/movie-track.mat';
load(filename)


trkdata = trk.data;
nanbouts = {};

for i=1:size(trkdata,1)
   currdata = trkdata(i,:,1);
   currdatanan = isnan(currdata);
   [starts,ends] = get_interval_ends(currdatanan);
   nanbouts.starts{i} = starts;
   nanbouts.ends{i} = ends;
   boutlengths = ends-starts;
   nanbouts.lengths{i} = boutlengths;
end

% write csv
[a,b] = fileparts(filename);
filesavename = fullfile(a,'nanbouts.csv');
% fid = fopen(filesavename,'w');

% fprintf(fid, '%s \n',filename);

for i = 1:numel(nanbouts.starts)
    for j = 1:numel(nanbouts.starts{i})
%     fprintf(fid, '%d, %d, %d, %d \n',i,nanbouts.starts{i}(j),nanbouts.ends{i}(j),nanbouts.lengths{i}(j));
    fprintf( '%d, %d, %d, %d \n',i,nanbouts.starts{i}(j),nanbouts.ends{i}(j),nanbouts.lengths{i}(j))

    end 

end 

% fclose(fid);

%% symbolic copy expdirs for testing registration
explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T211624_BLUE', ...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T211755_GRN', ...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T211859_RED', ...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T212543_RED-oval', ...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T212752_GRN-oval', ...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T212914_BLUE-oval'};

ignorefiles = {'ledregistrationimage.png','registrationimage.png','registered_trx.mat','registrationdata.mat','registrationdata.txt','indicatordata.mat'};
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers_try2';
for i = 1:numel(explist)
    expdir = explist{i};
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles)
end
%% symbolic copy expdirs for testing locomotion with new plate IDs
% 1b,1a,2b, 2a,4a,3b
 explist = {'/groups/branson/bransonlab/flydisco_data/pilot24nonLED_JHS_K_85321_RigA_20210313T074451', ...
     '/groups/branson/bransonlab/flydisco_data/pilot29nonLED_JHS_K_85321_RigA_20210312T094943', ...
     '/groups/branson/bransonlab/flydisco_data/pilot29nonLED_JHS_K_85321_RigB_20210315T145016', ...
     '/groups/branson/bransonlab/flydisco_data/pilot29nonLED_JHS_K_85321_RigB_20210312T093513', ...
     '/groups/branson/bransonlab/flydisco_data/pilot29nonLED_JHS_K_85321_RigD_20210312T095812', ...
     '/groups/branson/bransonlab/flydisco_data/pilot29nonLED_JHS_K_85321_RigD_20210315T140458'};
     
 ignorefiles = {'ledregistrationimage.png','registrationimage.png','registered_trx.mat','registrationdata.mat','registrationdata.txt','indicatordata.mat'};
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/testlocomotion_registration';
for i = 1:numel(explist)
    expdir = explist{i};
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles)
end    
     
%% symbolic copy expdirs for testing processing on Ryo's RGB experiments
explist = {'/groups/branson/bransonlab/flydisco_data/pilot29RGB_JHS_K_85321_RigB_20210315T163334', ...
        '/groups/branson/bransonlab/flydisco_data/pilot29RGB_EXT_VGLUT-GAL4_RigA_20210315T143243', ...
        '/groups/branson/bransonlab/flydisco_data/locomotionGtACR1_24_RGB_EXT_VGLUT-GAL4_RigA_20210226T095136'};
ignorefiles = {};    
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210319_RGBRyotestexperiments';
     
for i = 1:numel(explist)
    expdir = explist{i};
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles)
end    
    
%% grab a movie frame

% moviefile = '/groups/branson/bransonlab/alice/temp_bubbledata/LEDoffcenter/bias_video_v002.ufmf';
moviefile = '/groups/branson/bransonlab/alice/temp_bubbledata/LEDoffcenter/centered_carrot.ufmf';
moviefile = '/groups/branson/bransonlab/flydisco_data/locomotionGtACR1_24_RGB_EXT_VGLUT-GAL4_RigA_20210226T095136/movie.ufmf'

frame = 10;
readframe = get_readframe_fcn(moviefile);
frameimage = readframe(frame);
imshow(frameimage)

%% find screen_types to make settings dirs for Ryo's RGB experiments and which are already tracked

% explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/explist_LED_tobeprocessed_20210319.txt','%s');
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/explist_dextrosefood.txt','%s');
mdfile = 'Metadata.xml';


for i = 1:numel(explist)
    expdir = explist{i};
    metadata = ReadMetadataFile(fullfile(expdir, mdfile));
    fprintf('%s \n %s \n',expdir, metadata.screen_type);
    screen_type{i} = metadata.screen_type;
    if exist(fullfile(expdir,'movie-track.mat'),'file')
        fprintf('tracked \n')
    end
end


%% symbolic copy experiments to run APT 
% 
ignorefiles = {};    

% explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS29662_RigA_20210324T133352', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS33430_RigB_20210325T133338', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigC_20210324T161114', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_JHS_K_85321_RigB_20210325T152834'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654'};

% rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210427_testingFlyTrackerMods';
explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS52985_RigA_20210503T152107'};
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210506_testingMetadataMods';



for i = 1:numel(explist)
    expdir = explist{i};
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles)
end  
%% hard copy movielist
% explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS29662_RigA_20210324T133352', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS33430_RigB_20210325T133338', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_EXT_VGLUT-GAL4_RigC_20210324T161114', ...
%     '/groups/branson/bransonlab/flydisco_data/VNC_JHS_K_85321_RigB_20210325T152834'};
% rootoutdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210330_APTtracking';
explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS42702_RigC_20210505T153129', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS43522_RigC_20210412T131721', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS43660_RigA_20210414T135538', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS44231_RigA_20210419T141517', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS46301_RigC_20210512T131504', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS48619_RigA_20210426T141945', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS49172_RigA_20210504T145545', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS49220_RigB_20210421T143507', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS51893_RigA_20210506T150106', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS52985_RigA_20210503T152107', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS60240_RigB_20210420T141142', ...
'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS65692_RigC_20210427T135618', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS52094_RigA_20210602T144112', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS52799_RigD_20210602T145135', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS62706_RigA_20210525T134347', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS64175_RigB_20210518T145000', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS64214_RigA_20210520T151504', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS64218_RigA_20210527T162508', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS64234_RigC_20210513T154200', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS65710_RigB_20210513T142817', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS65765_RigA_20210525T140224', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS66932_RigA_20210525T153325', ...
'/groups/dickson/dicksonlab/flydisco_data/VNC_JRC_SS67323_RigA_20210513T145544'};

rootoutdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking';
for i = 1:numel(explist)
    expdir = explist{i};
    [~,expname] = fileparts(expdir);
    outdir = fullfile(rootoutdir,expname);
    cmd = sprintf('cp -R %s %s', expdir, outdir);
    system(cmd)
end

%% symbolic copy for testing FlyBowl data with new registration data
% 20210402
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654'};
% explist = {'/groups/branson/bransonlab/flydisco_data/VNC_YNA_K_162984_RigD_20210405T152704'}
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/nochr_TrpA71G01_Unknown_RigA_20201216T162938'};
% explist = {'/groups/branson/bransonlab/from_tier2/fly_bubble/bubble_data/nochr_TrpA65F12_Unknown_RigD_20201216T175902'};
% explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS25523_RigB_20210324T144555'};
% explist = {'/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS29892_RigD_20210323T134743'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654'};
% rootoutdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210419_wingtrackingfeatures_testing';
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigC_20210408T124720'};
% rootoutdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210420_testingCtraxresultsmovie';
rootoutdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210424_APTtracking';
% 
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/Untracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigA_20210408T130721', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/Untracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigB_20210408T130659', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/Untracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigC_20210408T124720', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/Untracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigD_20210408T130906'};
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/noLED/locomotionGtACR1_29_nonLED_JHS_K_85321_RigA_20210227T113908'};
rootindata = '/groups/branson/bransonlab/flydisco_data';
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/RyoRequests/TrackwithAPT.txt','%s');
for i = 1:numel(explist)
%     expdir = explist{i};
    expdir = fullfile(rootindata,explist{i})
    SymbolicCopyExperimentDirectory(expdir,rootoutdir)
end 


%% make short video clip - fly jumping over another fly. 
% moviefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/VNC_JRC_SS29892_RigD_20210323T134743';
% startframe = 1625
% endframe = 1635

%% frame rate
% 20210413
moviefile =  '/groups/branson/bransonlab/flydisco_data/VNC_JRC_SS25451_RigD_20210322T142656/movie.ufmf';
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);

%% testing command line project stufff
% with iObj = lbl open in APT
lObj.movieFilesUnMacroizeAll();
tblLabels = lObj.labelGetMFTableLabeled;
movs = lObj.movieFilesAllFull;
trxs = lObj.trxFilesAllFull;

%% look at nframes for VNC screen
% used to
% /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/Debug_getExperimentDirs_FlyBubbleRGB_Becca.m
% to grab all experiments and added nframes
%tempdate = datenum({expdirstruct.date},'yyyymmdd');
tempdate = str2double({expdirstruct.date});
h = figure;
ax =  plot(tempdate,[expdirstruct.nframes],'.');
ylim([50000,60000])

%% calculate dt cut offs to use in perframe stat calculations
moviefilelist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS42702_RigC_20210505T153129/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43522_RigC_20210412T131721/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43660_RigA_20210414T135538/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS44231_RigA_20210419T141517/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS46301_RigC_20210512T131504/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS48619_RigA_20210426T141945/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49172_RigA_20210504T145545/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49220_RigB_20210421T143507/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS51893_RigA_20210506T150106/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52094_RigA_20210602T144112/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52799_RigD_20210602T145135/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52985_RigA_20210503T152107/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS60240_RigB_20210420T141142/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS62706_RigA_20210525T134347/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64175_RigB_20210518T145000/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64214_RigA_20210520T151504/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64218_RigA_20210527T162508/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64234_RigC_20210513T154200/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65692_RigC_20210427T135618/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65710_RigB_20210513T142817/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65765_RigA_20210525T140224/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS66932_RigA_20210525T153325/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS67323_RigA_20210513T145544/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigA_20210519T153103/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigB_20210602T135902/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigC_20210604T152231/movie.ufmf',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigD_20210526T155055/movie.ufmf'};
for i = 1:numel(moviefilelist)
    i
    moviefile = moviefilelist{i};
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
dt = diff(headerinfo.timestamps);
maxdt(i)= max(dt);
mindt(i)= min(dt);
temp = sort(dt,'descend');
temp(1:5)
end
%% temporal registration
% loaded the seconds_fliesloaded from both metadata spreadsheets
% prctile(A,[50,75,90,93,95,99,99.5])
% ans =
%     1.7480    2.2217    3.2642    4.3869    6.7728   27.3700   29.9360
% sorted 
% 
%        NaN
%        NaN
%   458.2870
%   453.8990
%   331.3960
%   232.0640
%   153.6890
%    53.0720
%    37.3400
%    33.5900
% decided to have max cut off to fail
% trim movies to nframes param
%% pulled metadata from dickson and branson 
% combined expdirs
dickson = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_dicksonlab_20210805.mat');
branson = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_bransonlab_20210805.mat');

metadata = [branson.metadata,dickson.metadata];

save('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_ALL_20210805.mat','metadata');


%% make list for changing screen_type of 0212 led protocol experiments

savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirlist_0212ledprotocol.txt';
fid = fopen(savefile,'w');
for i = 1:numel(metadata_0212LEDprotocol)
    fprintf(fid, '%s\n', metadata_0212LEDprotocol(i).file_system_path);
end
fclose(fid);

%% pull timestamps for all movies 20210802
% load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_ALL_20210729.mat
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_ALL_20210805.mat
% for i = 1:numel(metadata)
%     expdir = metadata(i).file_system_path;
%     if exist(fullfile(expdir,'movie.ufmf'),'file')
% %         [~,~,fid,headerinfo] = get_readframe_fcn(fullfile(expdir,'movie.ufmf'));
%         headerinfo = ufmf_read_header(fullfile(expdir,'movie.ufmf'));
%         metadata(i).movielength = headerinfo.timestamps(end);
% %         fclose(fid);
%     else 
%         metadata(i).movielength = nan;
%     end
% end
[a,b,c] = unique({metadata.led_protocol});
metadata_0212LEDprotocol = metadata(c==1);
metadata_0331LEDprotocol = metadata(c==2);
save('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_ALL_20210805.mat','metadata','metadata_0212LEDprotocol','metadata_0331LEDprotocol');
% plotting 
figure, plot(datenum({metadata_0331LEDprotocol.date},'yyyymmdd'),[metadata_0331LEDprotocol.movielength],'.')
figure, plot(datenum({metadata_0212LEDprotocol.date},'yyyymmdd'),[metadata_0212LEDprotocol.movielength],'.')

%% find expdir with trx birth and death for testing FlyDiscoRegistration
expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS42702_RigC_20210505T153129', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43522_RigC_20210412T131721', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS43660_RigA_20210414T135538', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS44231_RigA_20210419T141517', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS46301_RigC_20210512T131504', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS48619_RigA_20210426T141945', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49172_RigA_20210504T145545', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS49220_RigB_20210421T143507', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS51893_RigA_20210506T150106', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52094_RigA_20210602T144112', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52799_RigD_20210602T145135', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS52985_RigA_20210503T152107', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS60240_RigB_20210420T141142', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS62706_RigA_20210525T134347', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64175_RigB_20210518T145000', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64214_RigA_20210520T151504', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64218_RigA_20210527T162508', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS64234_RigC_20210513T154200', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65692_RigC_20210427T135618', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65710_RigB_20210513T142817', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS65765_RigA_20210525T140224', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS66932_RigA_20210525T153325', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_JRC_SS67323_RigA_20210513T145544', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigA_20210519T153103', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigB_20210602T135902', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigC_20210604T152231', ...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210621_APTTracking/VNC_YNA_K_162984_RigD_20210526T155055'};

for i= 1:numel(expdirs)
    expdir = expdirs{i}
    load(fullfile(expdir,'registered_trx.mat'));
    numel(trx)
    [trx.nframes]
end
%% find unreadable metadata
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirlist_0212ledprotocol.txt','%s');
for i = 1:numel(explist)
    expdir = explist{i};
    metadatafile = fullfile(expdir,'Metadata.xml');
    fprintf('opening metadata %s:',expdir)
    try
        ReadMetadataFile(metadatafile);
        fprintf('success')
    catch
        fprintf('fail')
        
    end
    fprintf('\n')
end
%% should we auto fail dead or damaged videos?
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_bransonlab_20210805.mat
trajnum = [metadata.trajnum];
dead = [metadata.num_flies_dead];
damaged = [metadata.num_flies_damaged];

%%%%% change dead/damaged
[a,b,c] = unique(dead);

% plot
h = figure;
for i = 1:numel(a)
    subplot(1,numel(a),i)
    idx = [];
    idx = find(c == i);
    
    hist(trajnum(idx))
    trajnum(idx)
end
% ratios %%%%% change dead/damaged
for i = 1:numel(a)
    idx = [];
    idx = find(c == i);
    
    
    idx2 = [];
    idx2 = find(trajnum(idx) > 12);
    
    fprintf('%i dead flies: %04d : %04d experiments (<= 12 trajs: >12 trajs) \n',a(i), numel(idx)-numel(idx2),numel(idx2))
    
end
%% check for experiments that haven't run at all
load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_metadata_VNC_ALL_20210805.mat','metadata');
fid = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/experiments_without_tracker_data_and_notaborted.txt','w');
for i = 1:numel(metadata)
    expdir = metadata(i).file_system_path;
    trackerfilestr = 'movie-track.mat';
    trackerfile = fullfile(expdir,trackerfilestr);
    abort = fullfile(expdir,'ABORTED');
    if ~exist(trackerfile,'file') && ~exist(abort,'file')
        fprintf(fid,'%s\n',expdir);
    end
end

%% completed params
params = ReadParams('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings/20210531_flybubble_LED_ARtestingperframe/automatic_checks_complete_params.txt');
for i = 1:numel(params.required_files)
    fprintf('%s   %s \n', params.required_files{i},params.file_categories{i})
end

%% check on status of dickson lab reprocess including tracking
% retracked 8/9, but matlab crashed before caboose ran
% reran to do caboose step 8/16 but 5 movies started tracking then
% (autochecks incoming 'on' rather than 'force')
explist = textread('/groups/dickson/dicksonlab/Alice/RERUNposttracking_explist_dickson20210809T085057.txt','%s');
fid = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/metadata_dicksonRERUNtracking_20210816.csv','w');
autocheckscompletefilestr = 'automatic_checks_complete_results.txt';
trackingfilestr= 'movie-track.mat';
fprintf(fid, 'expdir, tracking file exist, tracking file date, ACC file exist, ACC date, ACC autopf, ACC category \n');
for i = 1:numel(explist)
    expdir = explist{i};
    fprintf(fid, '%s, ',expdir);
    trackfile = fullfile(expdir,trackingfilestr);
    if exist(trackfile,'file')
        tmp = [];
        tmp = dir(trackfile);
        filedate = datestr(tmp.date);
        fprintf(fid, '1, %s, ', filedate);
    else
        fprintf(fid, '0, Nan, ');
    end
        
    ACCfile = fullfile(expdir,autocheckscompletefilestr);
    if exist(ACCfile,'file')
        tmp2 =[];
        tmp2 = dir(ACCfile);
        ACCdate = datestr(tmp2.date);
        fprintf(fid,'1, %s, ',ACCdate);
        autocheckout = ReadParams(ACCfile);
        autopf = autocheckout.automated_pf;
        fprintf(fid,'%s, ',autopf);
        if autopf == 'F'
            fprintf(fid,'%s, \n',   autocheckout.automated_pf_category);
        else
            fprintf(fid,'\n');
        end
        
    else
        fprintf(fid, '0, Nan, \n');
        
    end
end
fclose(fid);
%% parse auto_pfs
 expdirstruct = metadata;
idx = strcmp({expdirstruct.automated_pf},'F');
temp2 = {expdirstruct(idx).automated_pf_category};
[d,e] = findgroups(temp2);
counts = histcounts(d);
for i = 1:numel(e)
fprintf('%2d, %s \n', counts(i),e{i})
end
%% count line expdirs
temp2 = cellstr({expdirstruct.line});
[d,e] = findgroups(temp2);
counts = histcounts(d);
for i = 1:numel(e)
fprintf('%2d, %s \n', counts(i),e{i})
end
%% combine VNC experiment metadata
dlab = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_autochecks_VNC_dicksonlab_20210818.mat');
blab = load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_autochecks_VNC_bransonlab_20210818.mat');

metadata = [blab.metadata,dlab.metadata];
save('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_autochecks_VNC_combined_20210818.mat','metadata','blab','dlab')

expdirtable = struct2table(metadata);
writetable(expdirtable,'/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_autochecks_VNC_combined_20210818.csv','Delimiter','tab');

%% find non-VNC expdirs
% pulled all data dirs in bransonlab with ~FlyDiscoAnalysis/Debug_getExperimentDirs_FlyBubbleRGB.m
% expdirstruct(1).file_system_path ans =    '.DS_Store'

metadata = expdirstruct(2:end);
 
scrntype = cellstr({metadata.screen_type});

unique(scrntype)'

%     {'non_olympiad_dickson_VNC'                       }
%     {'non_olympiad_dickson_led5secVNC'                }
%     {'non_olympiad_dickson_locomotionGtACR1'          }
%     {'non_olympiad_dickson_locomotionGtACR1_24_RGB'   }
%     {'non_olympiad_dickson_locomotionGtACR1_24_nonLED'}
%     {'non_olympiad_dickson_locomotionGtACR1_29_RGB'   }
%     {'non_olympiad_dickson_locomotionGtACR1_29_nonLED'}
%     {'non_olympiad_dickson_pilot24RGB'                }
%     {'non_olympiad_dickson_pilot24nonLED'             }
%     {'non_olympiad_dickson_pilot29RGB'                }
%     {'non_olympiad_dickson_pilot29nonLED'             }
%     {'non_olympiad_dickson_pilotPTR'                  }
%     {'non_olympiad_dickson_pilotRGB16'                }

othertypes = {'non_olympiad_dickson_locomotionGtACR1','non_olympiad_dickson_locomotionGtACR1_24_RGB','non_olympiad_dickson_locomotionGtACR1_24_nonLED', ...
    'non_olympiad_dickson_locomotionGtACR1_29_RGB','non_olympiad_dickson_locomotionGtACR1_29_nonLED','non_olympiad_dickson_pilot24RGB', ...
    'non_olympiad_dickson_pilot24nonLED','non_olympiad_dickson_pilot29RGB','non_olympiad_dickson_pilot29nonLED','non_olympiad_dickson_pilotPTR', ...
    'non_olympiad_dickson_pilotRGB16'};
idx1 = strcmp(scrntype,'non_olympiad_dickson_VNC');
idx2 = strcmp(scrntype,'non_olympiad_dickson_led5secVNC');
idx = idx1 + idx2;

metadata(~idx).file_system_path;
explist = cellstr({metadata(~idx).file_system_path});

fid = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/nonscreen_expdirs_20210913.txt','w');
fprintf(fid,'%s, \n',metadata(~idx).file_system_path);
fclose(fid);
save('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/nonscreen_expdirs_20210913','explist');
%% check for aborted with flag but not file. 
idxa = find([metadata.flag_aborted]== 1);
% metadata(idxa).flag_aborted;
metadata(idxa).file_system_path;
explist = cellstr({metadata(idxa).file_system_path});
fid = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/aborted_flag_expdirs_20210913.txt','w');
fprintf(fid,'%s \n',metadata(idxa).file_system_path);
fclose(fid);
save('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/aborted_flag_expdirs_20210913','explist');
%% convert Ryo's tracking request list 11/4/2021
datadir = '/groups/branson/bransonlab/flydisco_data';
explistog = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/expdirlist_toAPTtrack.txt','%s');
fid = fopen('/groups/branson/home/robiea/Projects_data/FlyDisco/expdirlist_toAPTtrack_20211104_linuxpathsFIXED.txt','w');
for i = 1:numel(explistog)
    currdirpath = fullfile(datadir,explistog{i});
    fprintf(fid,'%s\n',currdirpath);
    
    
end
fclose(fid);
%% symbolic copy expdirs for testing processing on Ryo's RGB experiments
cd /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/Katies_expdirlist_20211112.txt','%s');
ignorefiles = {};    
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211111_testingmetadatachangesForKatie';
     
for i = 1:numel(explist)
    expdir = explist{i};
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles)
end  

%% run JAABADetect from /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe
cd /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/Katie_expditlist_20211113_localpath.txt','%s');
auto_checks = 'automatic_checks_complete_results.txt';
jabfiles = {'/groups/branson/home/bransonlab/JAABAProjectFiles/FlyBubble/FeAggBubDisco8.jab'};
for i = 1:numel(explist)
    expdir = explist{i}
%     accfile = fullfile(expdir,auto_checks);
%     if exist(accfile,'file')
%         AutoChecksComplete = ReadParams(accfile)
%         if strcmp(AutoChecksComplete.automated_pf,'P');
        JAABADetect(expdir,'jabfiles',jabfiles,'forcecompute',false);

%         end
%     end
    
    
end
%% how many partial trajectories
explist = textread('/groups/branson/home/robiea/Projects_data/Katie/NorpA_Bubble/bubble_explist20211117_NorpAaIPg.txt','%s');
for i = 1:numel(explist)
    expdir = explist{i}
    load(fullfile(expdir,'registered_trx.mat'))
    trx.off
    
end
%% checking for jitter in rigs 3/23/2022
explist = {'/groups/branson/bransonlab/flydisco_data/NorpA_JHS_K_85321_RigB_20210903T074659',...
    '/groups/branson/bransonlab/flydisco_data/DecreAggTen2_JRC_SS36564_RigB_20220215T082557',...
'/groups/branson/bransonlab/flydisco_data/DecreAggTen2_JRC_SS36564_RigB_20220215T080215',...
'/groups/branson/bransonlab/flydisco_data/DecreAggTen2_JRC_SS36564_RigA_20220215T085842',...
'/groups/branson/bransonlab/flydisco_data/DecreAggR53_JRC_SS36564_RigC_20220215T085138',...
'/groups/branson/bransonlab/flydisco_data/DecreAggR53_JRC_SS36564_RigC_20220215T082833',...
'/groups/branson/bransonlab/flydisco_data/DecreAggR53_JRC_SS36564_RigA_20220215T080013',...
'/groups/branson/bransonlab/flydisco_data/DecreAgg53_JRC_SS36564_RigC_20220215T074958',...
'/groups/branson/bransonlab/flydisco_data/DecreAgg53_JRC_SS36564_RigB_20220215T084926',...
'/groups/branson/bransonlab/flydisco_data/DecreAgg53_JRC_SS36564_RigA_20220215T083529',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigC_20220215T083935',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigC_20220215T081640',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigB_20220215T081413',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigA_20220215T074630',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_pBDPLexAp65U_attP40_RigB_20220215T083728',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_pBDPLexAp65U_attP40_RigB_20220215T074748',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigA_20220215T081222',...
'/groups/branson/bransonlab/flydisco_data/CsChr53_JRC_SS36564_RigA_20220215T084712',...
'/groups/branson/bransonlab/flydisco_data/DecreAggTenR2_JRC_SS36564_RigC_20220215T080423',...
'/groups/branson/bransonlab/flydisco_data/DecreAggTenR2_JRC_SS36564_RigB_20220215T090045',...
'/groups/branson/bransonlab/flydisco_data/DecreAggTenR2_JRC_SS36564_RigA_20220215T082353',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220217_flybubble_TrpA/TrpAFemale2_GMR_Eb5_vk5_RigB_20220216T085536',...
'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220217_flybubble_TrpA/TrpAMale2_GMR_71G01_JK73A_RigB_20220216T074837'};
% explist = {'/groups/branson/bransonlab/flydisco_data/NorpA_JHS_K_85321_RigB_20210903T074659'};

  analysis_protocol = '20210806_flybubble_LED_NorpA';
  figure
  Nexp = numel(explist);

  for i = 1:Nexp
  expdir = explist{i};
  [dx,dy,xs,ys] = MeasureJitter(expdir,'analysis_protocol',analysis_protocol,'endoff',300);
  idx = find(~isnan(xs));
  meanx = mean(xs(idx));
  meany = mean(ys(idx));
  dmx = meanx - xs(idx);
  dmy = meany - ys(idx);
  subplot(12,2,i)
  [~,expname] = fileparts(expdir);
  title(expname,'interpreter','none');
  hold on
  plot(dmx,'r')
  plot(dmy,'b')
  
  end
  
  %plot

  figure
  plot(xs(idx))
  
%% explore metadata 
% load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20220406.mat
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNCall_20221102.mat
% /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC_20220406.csv

% restrict to VNC and VNC2 ( also contains non_olympiad_dickson_led5secVNC)
expdirstruct = metadata;
idx1 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC');
idx2 = strcmp({expdirstruct.screen_type},'non_olympiad_dickson_VNC2');
idxall = idx1 + idx2;
expdirstruct = expdirstruct(logical(idxall));

[a,b,c] = unique({metadata.automated_pf});
passnum = numel(find(c == 2)); %3155
failnum = numel(find(c==1));%182
% nannum = numel(find(c==2));

metadata_fail = metadata(c==1);
[aa,bb,cc] = unique({metadata_fail.automated_pf_category})



[d, id] = findgroups({metadata_fail.automated_pf_category});
counts = histcounts(d)
%% line counts 20221103
% [a,b,c] = unique({metadata.line});
[a,b,c] = unique({metadata.automated_pf});
metadata_pass = (metadata( c==2));
[d, id] = findgroups({metadata_pass.line});
counts = histcounts(d,numel(id))

%% reprocessing for Katie
load /groups/branson/home/robiea/Projects_data/Katie/TrpAforReprocessing.mat

explist = {metadata.file_system_path};

toprocessdir =  '/groups/branson/bransonlab/flydisco_data/to-process';
%ln -s source_file_or_directory_name  softlink_name
for i = 1:numel(explist)
    expdir = explist{i};
    [~,expname] = fileparts(expdir);
    outdir = fullfile(toprocessdir,expname);
    cmd = sprintf('ln -s %s %s', expdir, outdir);
    system(cmd)
end

%% how did ryo organize 3 days?
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNCOCt2021_3perweek_20220428.mat

% [c,ia,ic] = unique({metadata.line});
dates = {metadata.date};
lines = {metadata.line};
dates_num = datenum(dates,'yyyymmdd');
idx = find(dates_num < datenum('20211010','yyyymmdd'));
% figure, plot(dates_num(idx),'.')

for i = 1:numel(idx)
sprintf('%s,%s\n',lines{idx(i)},datestr(dates_num(idx(i)),'yyyymmdd'))
end

%% check protocol vs detected 
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220503_flybubble_VNC2_testLEDcountlength/VNC2_YNA_K_162984_RigB_20220427T115238';
% expdir = '/groups/branson/bransonlab/alice/temp_data/disco_LEDtesting/VNC2_JRC_SS90575_GtACR1_RigD_20220506T173925';
% expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220507T182908';
% expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_EXT_VGLUT-GAL4_RigA_20220503T095111';

load(fullfile(expdir,'protocol.mat'));
load(fullfile(expdir,'indicatordata.mat'));
moviefile = (fullfile(expdir,'movie.ufmf'));
[~,~,~,headerinfo] = get_readframe_fcn(moviefile);
% turn 3 protocol into 1 color protocol
RGBprotocol = protocol;
protocol = downmixProtocolIfNeeded(metadata, RGBprotocol) ;

% check stim count
stimcount_expected  = sum(protocol.iteration);
stimcount_detection = max(numel(indicatorLED.startframe),numel(indicatorLED.endframe));
if ~stimcount_expected == stimcount_detection
    fprintf('Incorrect number of bouts\n')
else
    fprintf('Correct number of bouts\n')
end


% check stim duration in time
% find expected stim durations
stimdurs_expected = [];
for i = 1:max(protocol.stepNum)
   currstimMS = protocol.pulseNum(i)*protocol.pulsePeriodSP(i) - (protocol.pulsePeriodSP(i) - protocol.pulseWidthSP(i)); 
   currstep_stimdurs_expected = repmat(currstimMS,1,protocol.iteration(i));
   stimdurs_expected = [stimdurs_expected,currstep_stimdurs_expected];
end

%assuming LEDs off at beggining and end of experiment

stimcnt = 0;
stimdiffMS = nan(1,sum(protocol.iteration));
stimdiffFrm = nan(1,sum(protocol.iteration));
for i = 1:max(protocol.stepNum)
   for  j = 1:protocol.iteration(i)
    stimcnt = stimcnt + 1;
    currstimdur = [];
    currstimdur = (indicatorLED.endtimes(stimcnt)-indicatorLED.starttimes(stimcnt))*1000;
    currstimdurFrm = (indicatorLED.endframe(stimcnt)-indicatorLED.startframe(stimcnt));
    stimdiffMS(stimcnt) = stimdurs_expected(stimcnt)-currstimdur;
    stimdiffFrm(stimcnt) = stimdurs_expected(stimcnt)/1000*150 - currstimdurFrm;
   end
    
end

% set threshold of frames to be off
frmnum4threshold = 1;
frmlengthMS = 1/150*1000;
stimdurthreshold = frmnum4threshold * frmlengthMS;

if any(find(abs(stimdiffMS) > stimdurthreshold ))
    fprintf('Bout length incorrect for stimulus:,%s \n', num2str(find(abs(stimdiffMS) > stimdurthreshold )))
else 
    fprintf('Correct bout lengths\n')
    
end

%% compare protocols
cd '/groups/branson/home/robiea/Code_versioned/FlyBowlDataCapture'
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220503_flybubble_VNC2_testLEDcountlength/VNC2_pBDPGAL4U_RigD_20220507T182908';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220503_flybubble_VNC2_testLEDcountlength/VNC2_YNA_K_162984_RigB_20220427T115238';
% expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220508T180927';
% expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220508T183056';
% expdir = '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220508T182327';
explist = {'/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigB_20220509T143008', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigB_20220509T141645', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigA_20220509T143008', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigA_20220509T141411', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220509T143013', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigD_20220509T141756', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigC_20220509T143013', ...
    '/groups/branson/bransonlab/flydisco_data/VNC2_pBDPGAL4U_RigC_20220509T141539'};

for i = 1:numel(explist)
    expdir = explist{i};
    load(fullfile(expdir, 'protocol.mat'));
    load(fullfile(expdir, 'indicatordata.mat'));
    moviefile = (fullfile(expdir,'movie.ufmf'));
    [~,~,~,headerinfo] = get_readframe_fcn(moviefile);
    timestamps = headerinfo.timestamps;
    
    hfig = figure;
    hax = gca;
    %from FBDC code base
    DisplayStimulusProtocol(protocol,'hax',hax)
    hold on
    plot(hax,timestamps,indicatorLED.indicatordigital)
    xlim([0,max(timestamps)])
end
%% copy and delete expdirs
explist = {'VNC2_pBDPGAL4U_RigD_20220507T182908',...
'VNC2_pBDPGAL4U_GtACR1_RigD_20220507T184141',...
'VNC2_pBDPGAL4U_GtACR1_RigB_20220507T185552',...
'VNC2_pBDPGAL4U_GtACR1_RigB_20220507T190633',...
'VNC2_pBDPGAL4U_GtACR1_RigD_20220507T191701',...
'VNC2_pBDPGAL4U_GtACR1_RigD_20220507T192407',...
'VNC2_pBDPGAL4U_RigD_20220508T180927',...
'VNC2_pBDPGAL4U_RigD_20220508T182327',...
'VNC2_pBDPGAL4U_RigD_20220508T183056',...
'VNC2_pBDPGAL4U_RigA_20220509T141411',...
'VNC2_pBDPGAL4U_RigC_20220509T141539',...
'VNC2_pBDPGAL4U_RigB_20220509T141645',...
'VNC2_pBDPGAL4U_RigD_20220509T141756',...
'VNC2_pBDPGAL4U_RigA_20220509T143008',...
'VNC2_pBDPGAL4U_RigB_20220509T143008',...
'VNC2_pBDPGAL4U_RigC_20220509T143013',...
'VNC2_pBDPGAL4U_RigD_20220509T143013'};

ignorefiles = {};    
rootoutputdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220610_testingFBDCfix';
%copy off dm11     
for i = 1:numel(explist)
    expdir = fullfile('/groups/branson/bransonlab/flydisco_data',explist{i});
    SymbolicCopyExperimentDirectory(expdir,rootoutputdir,'ignorefiles',ignorefiles,'dosoftlink',false)
end 
%delete off flydisco_data (have to run as bransonlab
for i = 1:numel(explist)
    expdir = fullfile('/groups/branson/bransonlab/flydisco_data',explist{i});
    cmd= sprintf('rm -r %s\n',expdir)
    system(cmd)
end
%% find VNC2 expdir < 20220523 to remake ctrax results videos   
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_VNC2_20220611.mat
date = {metadata.date};
for i = 1:numel(date)
idx(i) = str2num(date{i})< 20220523;
end
% sort(date(idx))
% sort(date(~idx))
explist = {metadata(idx).file_system_path}';
% explist = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220608_testpipelinefiledeletion/VNC2_JRC_SS83407_RigB_20220419T092554'};
% remove ctrax results movie and make softlink in to-process folder
todeletewildcard = 'ctrax_results_movie_*.mp4';


toprocessdir =  '/groups/branson/bransonlab/flydisco_data/to-process';
%ln -s source_file_or_directory_name  softlink_name
for i = 1:numel(explist)
    [~,expdir] = fileparts(explist{i});
    filestr = fullfile(explist{i},todeletewildcard);
    if ~isempty(dir(filestr))
        resultsmoviename = dir(filestr);
                delete(fullfile(explist{i},resultsmoviename.name))
                if ~isempty(dir(fullfile(explist{i},resultsmoviename.name)))
                    fprintf('not deleted %s\n',resultsmoviename.name)
                    return
                end
    end

    [~,expname] = fileparts(expdir);
    outdir = fullfile(toprocessdir,expname);

    cmd = sprintf('ln -sfn %s %s', explist{i}, outdir);

%     cmd = sprintf('ln -s %s %s', expdir, outdir);
    system(cmd)
end
%% add manual_failure to VNC expdir data pull
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116.mat
savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116_addedManualFail';
for i= 1:numel(expdirstruct2)
    failfile = fullfile(expdirstruct2(i).file_system_path,'manual_fail.txt');
    if exist(failfile,'file')
        expdirstruct2(i).manual_fail = 'F';
        manual_fail_category = textread(failfile,'%s','delimiter','\n');
        expdirstruct2(i).manual_fail_category = manual_fail_category;
%         fprintf('%s\n',failfile)
    end
end
save(savefile,'expdirstruct2');
%% print out tsv
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116_addedManualFail.mat
savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116_addedManualFail';

expdirstruct = expdirstruct2;
fid = fopen([savefile,'.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t  %s\t  %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','trajnum','notes_tech','notes_behav','automated_pf','auto_pf_category','manual_fail','manual_fail_category');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  trajnum = expdirstruct(i).trajnum;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;
  autopfcat = expdirstruct(i).automated_pf_category;
  manf = expdirstruct(i).manual_fail;
  manfcat = cell2str(expdirstruct(i).manual_fail_category);


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t  %d\t  %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,trajnum,notestech,notesbeh,autopf,autopfcat,manf,manfcat);
end

fclose(fid);
%% look at tracking for trajectory number 12+ 
% all VNC and VNC2 data through 2022
% explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/hightrajnum_explist.txt','%s');

% VNC2 data from 2023 
explist = textread('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirlist_hightrajnum_2023.txt','%s');
rootdatadir = '/groups/branson/bransonlab/flydisco_data';

for i = 1:numel(explist)
    if exist(fullfile(rootdatadir,explist{i},'manual_fail.txt'),'file')
        i
        disp('manual fail')
        continue
    elseif exist(fullfile(rootdatadir,explist{i},'automatic_checks_complete_results.txt'),'file')
        autochecks = ReadParams(fullfile(rootdatadir,explist{i},'automatic_checks_complete_results.txt'));
        if autochecks.automated_pf == 'F'||strcmp(autochecks.automated_pf,'NaN')
            i
            disp('auto fail')
            continue
        elseif isempty(dir(fullfile(rootdatadir,explist{i},'*.mp4')))
            i
            disp(explist{i})
            playfmf('filename',fullfile(rootdatadir,explist{i},'movie.ufmf'))
            uiwait
        else
            i
            moviename = dir(fullfile(rootdatadir,explist{i},'ctrax*.mp4'));
            disp(fullfile(rootdatadir,explist{i},moviename.name))
            implay(fullfile(rootdatadir,explist{i},moviename.name))
            uiwait
        end
    elseif isempty(dir(fullfile(rootdatadir,explist{i},'ctrax*.mp4')))
        i
        disp(explist{i})
        playfmf('filename',fullfile(rootdatadir,explist{i},'movie.ufmf'))
        uiwait
    else
        i
        moviename = dir(fullfile(rootdatadir,explist{i},'*.mp4'));
        disp(fullfile(rootdatadir,explist{i},moviename.name))
        implay(fullfile(rootdatadir,explist{i},moviename.name))
        uiwait
    end

end

%% compile numbersfor bad fly numbers by experimenter
load /home/robiea@hhmi.org/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_allflydisco_20230301_allVNC.mat
% problem - empty categories are numbers
idx = strcmp({metadata.automated_pf},'F');
failed = {metadata(idx).automated_pf_category};

% two ways to count 
% categories
c = categorical(failed);
categories(c);
countcats(c);

% unique
[a,~,d] = unique(failed);
counts = histcounts(d,1:numel(a));

% experimenters of failed experiments
experimentersF = {metadata(idx).experimenter};
c = categorical(experimenters);
categories(c);
countcats(c);
% total experiments 
experimenters = {metadata.experimenter};
c = categorical(experimenters);
categories(c)
countcats(c)'

% manual fails

idx2 = strcmp({metadata.manual_fail},'F');
mfailed = [metadata(idx2).manual_fail_category];
e = categorical(mfailed);
[aa,~,cc]= unique(e);
counts = histcounts(cc)
experimenters = {metadata(idx2).experimenter};
c = categorical(experimenters);
categories(c);
countcats(c);

%% dm11 flydisco_data clean up 9/1/2023 - explore data and make list for non_olympiad_dickson_led5secVNC
% pulled all directories in flydisco with
% /groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/Debug_getExperimentDirs_FlyBubbleRGB.m 

load('/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_allflydisco_20230831_20230901T132911.mat')

% make a list of screen_type non_olympiad_branson_AmpRec - move to
% nearline

% make a list of screen_type 
[a,b,c] = unique(cellfun(@num2str,{expdirstruct.screen_type},'uni',0));
[counts, groupnames] = groupcounts(c)
for i = 1:numel(counts),fprintf('%s has %d expdirs \n',expdirstruct(b(groupnames(i))).screen_type,counts(i)),end

% list dates of VNC experiments
VNCidx = strcmp({expdirstruct.screen_type}, 'non_olympiad_dickson_VNC');
VNCexpdirstruct = expdirstruct(VNCidx);
[a,b,c] = unique(cellfun(@num2str,{VNCexpdirstruct.date},'uni',0));
a'
% find 20210322 in VNC
wrongscreen_typeIdx = strcmp({VNCexpdirstruct.date}, '20210322');
VNCexpdirstruct(wrongscreen_typeIdx).file_system_path


% list dates of VNC2 experiments
VNC2idx = strcmp({expdirstruct.screen_type}, 'non_olympiad_dickson_VNC2');
VNC2expdirstruct = expdirstruct(VNC2idx);
[a,b,c] = unique(cellfun(@num2str,{VNC2expdirstruct.date},'uni',0));
a'
% list dates of led5secVNC
VNCLED5idx = strcmp({expdirstruct.screen_type}, 'non_olympiad_dickson_led5secVNC');
VNCLED5expdirstruct = expdirstruct(VNCLED5idx);
[a,b,c] = unique(cellfun(@num2str,{VNCLED5expdirstruct.date},'uni',0));
a'
% make list of led5secVNC experiments to move to nearline 9/1/2023
filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/led5secVNC_cleanuplist.txt'
fid = fopen(filename,'w')
for i = 1:numel(VNCLED5expdirstruct)
    fprintf(fid,'%s\n',VNCLED5expdirstruct(i).file_system_path)
end
fclose(fid)

%% dm11 flydisco_data clean up 9/1/2023 and make a list of make a list of automated_pf = F or manual-fail = F for VNC and VNC2 screen_types
% for data from march 2021 - oct 2022
clear all
% load data for VNC before 2023 data added
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116_addedManualFail.mat
% restrict to VNC and VNC2
idx1 = strcmp({expdirstruct2.screen_type},'non_olympiad_dickson_VNC');
idx2 = strcmp({expdirstruct2.screen_type},'non_olympiad_dickson_VNC2');
idxall = idx1 +idx2;
expdirstruct2 = expdirstruct2(logical(idxall));
% % explore data
% [a,b,c] = unique(cellfun(@num2str,{expdirstruct2.automated_pf},'uni',0));
% counts = groupcounts(c)
% % f - 360
% % NaN - 379
% % p - 5121
% 
% % what are the NaNs
% NaNsIdx = strcmp({expdirstruct2.automated_pf},'NaN');
% NaNexpdirstruct = expdirstruct2(NaNsIdx);
% [a,b,c] = unique(cellfun(@num2str,{NaNexpdirstruct.manual_fail},'uni',0));
% counts = groupcounts(c)
% % none - 191
% % manual_failed - 188
% % from spreadsheet /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20221116_addedManualFail.tsv
% % not manually failed with auto_pf = NaN -> 185 are led5secVNC, 6 are
% % pipeline folders like to-be-processed

% make a list of automated_pf = F or manual-fail = F

idx3 = strcmp({expdirstruct2.automated_pf},'F');
autofailedexpdirstruct = expdirstruct2(idx3);
idx4 = strcmp({expdirstruct2.manual_fail},'F');
manualfailedexpdirstruct = expdirstruct2(idx4);

AFailList = {autofailedexpdirstruct.file_system_path};
MFailList = {manualfailedexpdirstruct.file_system_path};

AandMFaillist = [AFailList,MFailList];

filename = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/AandMFail_VNCandVNC2_cleanuplist.txt';
fid = fopen(filename,'w');
for i = 1:numel(AandMFaillist)
    fprintf(fid,'%s\n',AandMFaillist{i});
end
fclose(fid);

% manually removed
% /groups/branson/bransonlab/flydisco_data/TrpAFemale4_EXT_VGLUT-GAL4_RigD_20220713T094246
% (wrong screen_type) 

%% high trajectory number - 10/2023
% print out tsv from data pull with auto and manual fails (from logbook so
% far) 
load /groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/expdirs_2023_pulled202310032023_20231004T103003.mat
savefile = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline/allVNC_20231004_addedManualFail';

fid = fopen([savefile,'.tsv'],'w');

fprintf(fid,'%s\t %s\t %s\t %s\t  %s\t  %s\t %s\t %s\t %s\t %s\t %s \n','expname','date','datetime','linename','trajnum','notes_tech','notes_behav','automated_pf','auto_pf_category','manual_fail','manual_fail_category');

for i = 1:numel(expdirstruct)  
  [~,expname] = fileparts(expdirstruct(i).file_system_path);
  datestr = expdirstruct(i).date;
  datetimestr = expdirstruct(i).exp_datetime;
  linename = expdirstruct(i).line;
  trajnum = expdirstruct(i).trajnum;
  notestech = expdirstruct(i).notes_technical;
  notesbeh = expdirstruct(i).notes_behavioral;
  autopf = expdirstruct(i).automated_pf;
  autopfcat = expdirstruct(i).automated_pf_category;
  manf = expdirstruct(i).manual_fail;
  manfcat = expdirstruct(i).manual_fail_category; % code above added manual fail posthoc not from pull


%   notes = experiments_eddison(i).notes_curation;
  fprintf(fid,'%s\t %s\t %s\t %s\t  %d\t  %s\t %s\t %s\t %s\t %s\t %s \n',expname, datestr,datetimestr,linename,trajnum,notestech,notesbeh,autopf,autopfcat,manf,manfcat);
end

fclose(fid);

