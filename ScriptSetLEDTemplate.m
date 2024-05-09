%% set up paths
    modpath
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
% addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
% settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings/';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings/';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';


%% parameters
% FlyBubbeRGB
% analysis_protocol = '20210219_flybubble_dickson_RGBGtACR1testing';
% analysis_protocol = '20210311_flybuthanksbble_dickson_RGBGtACR1';
% expdirs = {'/groups/branson/bransonlab/alice/temp_bubbledata/singlecolormarkers/SmallOval/Green/pilot24nonLED_JHS_K_85321_GtACR1_RigC_20210317T214041'};

% FlyBowlRGB
% analysis_protocol = '20210329_flybubble_flybowloptoKatie_mingrig_flytracker'; % changed name to 20210329_flybubble_FlyBowlRGB_LED
% analysis_protocol = '20210329_flybubble_FlyBowlRGB_LED';
% analysis_protocol = '20220622_flybubbleRed_MBL';
% analysis_protocol = '20240124_multibubble_firsttry';
% analysis_protocol = '20240124_multibubble_secondtry';
analysis_protocol = '20240507_flybubble_LED_VNC3'

% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/KatieTestData/FlyBowlDisco_RGBonly_318/20210318T135921_rig1_flyBowl2__SS36564_CsChrim_KS_redonly_protocolRGB_0315_2'}
% adjusted camera height
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/KatieTestData/FlyBowlDisco_RGBonly_401/20210401T132850_rig1_flyBowl2__aIPgSS1UASCsChrimson_KS_redonly_protocolRGB_0315_2'};
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220617_testingFlyBubbleRed/CsChrSocial3_P1a_Unknown_RigF_20220616T153431'};
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20220622_testingFlyBubbleRedMBLhardware/CsChrSocial3_TK_RigF_20220622T173359'};
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/20240122T154843_rig2_flyBowl1_1xLwt_attp40_4stop1_3xDSCPLwt_attp40_3stop1_SlowRamp';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240214T095625_rig1_multibubble_CSx71G01_1_1_UAS_Chrimson_Venus_X_0070_SlowRamp';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240507_testingVNC3/VNC3_YNA_K_162984_RigD_20240507T182533';

datalocparamsfilestr = 'dataloc_params.txt';
dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));


%%
moviefilename = fullfile(expdir,dataloc_params.moviefilestr);

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
im = readfcn(1);
% for j = 1:100:(headerinfo.nframes/5)
for j = 1:10:(headerinfo.nframes)
    tmp = readfcn(j);
    idx = tmp > im;
    im(idx) = tmp(idx);
    if (mod(j,1000) == 0);
        disp(round((j/headerinfo.nframes)*100))
    end
end


%% grab a rectangle around the marker to make a template
mu = im;
figure(1);
clf;
imagesc(mu);
axis image;
r = getrect;
xlim = round(r(1)+[0,r(3)]);
ylim = round(r(2)+[0,r(4)]);
hold on;
plot(xlim([1,1,2,2,1])+.5*[-1,-1,1,1,-1],ylim([1,2,2,1,1])+.5*[-1,1,1,-1,-1],'w-','linewidth',2);

template = mu(ylim(1):ylim(2),xlim(1):xlim(2));
figure(2);
imagesc(template);
axis image;

%%%%%%%%% change output file name; copy into analysis protocol directory and edit registration_params.LEDMarkerType  %%%%%%%%%%%%%%%%
imwrite(uint8(template),'LEDTemplates_multiBubble.png');

%% distance from the corner 
% informative for setting registration_params.maxDistCornerFrac_LEDLabel (restricts template
% matching to smaller area)
x = mean(xlim);
y = mean(ylim);
nc = headerinfo.nc;
nr = headerinfo.nr;
corners = [1,1,nc,nc;1,nr,nr,1];
dcorner = sqrt(min((corners(1,:)-x).^2 + (corners(2,:)-y).^2));
imr = min(nr,nc);

dcornerfrac = dcorner / imr

% maxDistCornerFrac_BowlLabel = .12;