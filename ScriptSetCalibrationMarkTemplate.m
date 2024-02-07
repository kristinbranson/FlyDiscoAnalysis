%% set up paths
modpath
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
% addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
% settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings/';
% settingsdir =
% '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings/';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';

%% parameters

% % analysis_protocol = '20150717_flybubble_flybowlMing';
analysis_protocol = '20240124_multibubble_firsttry';
% expfile = '/groups/branson/home/robiea/Projects_data/Ming/explist_forregistration';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/20240122T154843_rig2_flyBowl1_1xLwt_attp40_4stop1_3xDSCPLwt_attp40_3stop1_SlowRamp';
datalocparamsfilestr = 'dataloc_params.txt';
dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));



%%  
moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
annfilename = fullfile(expdir,dataloc_params.annfilestr);
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
im = readfcn(1);

mu = im;

%% grab a rectangle around the marker to make a template

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

%%%%%%%%% change output file name; copy into analysis protocol directory and edit registration_params.bowlMarkerType  %%%%%%%%%%%%%%%%
imwrite(uint8(template),'CirleTemplate_FlyBowlMing.png');

%% distance from the corner
% informative for setting registration_params.maxDistCornerFrac_BowlLabel (restricts template
% matching to smaller area)
x = mean(xlim);
y = mean(ylim);
corners = [1,1,nc,nc;1,nr,nr,1];
dcorner = sqrt(min((corners(1,:)-x).^2 + (corners(2,:)-y).^2));
imr = min(nr,nc);
dcornerfrac = dcorner / imr;

% maxDistCornerFrac_BowlLabel = .12;