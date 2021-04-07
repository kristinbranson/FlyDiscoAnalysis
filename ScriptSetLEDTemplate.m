%% set up paths
    
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
% settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings/';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings/';


%% parameters
% FlyBubbeRGB
% analysis_protocol = '20210219_flybubble_dickson_RGBGtACR1testing';
% analysis_protocol = '20210311_flybuthanksbble_dickson_RGBGtACR1';
% expdirs = {'/groups/branson/bransonlab/alice/temp_bubbledata/singlecolormarkers/SmallOval/Green/pilot24nonLED_JHS_K_85321_GtACR1_RigC_20210317T214041'};

% FlyBowlRGB
% analysis_protocol = '20210329_flybubble_flybowloptoKatie_mingrig_flytracker'; % changed name to 20210329_flybubble_FlyBowlRGB_LED
analysis_protocol = '20210329_flybubble_FlyBowlRGB_LED';
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/KatieTestData/FlyBowlDisco_RGBonly_318/20210318T135921_rig1_flyBowl2__SS36564_CsChrim_KS_redonly_protocolRGB_0315_2'}
% adjusted camera height
expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/KatieTestData/FlyBowlDisco_RGBonly_401/20210401T132850_rig1_flyBowl2__aIPgSS1UASCsChrimson_KS_redonly_protocolRGB_0315_2'};


datalocparamsfilestr = 'dataloc_params.txt';
dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));

% expdirs = importdata(expfile);
% expdirs(cellfun(@isempty,expdirs)) = [];
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/FlyBubbleRGB/locomotionGtACR1_24_RGB_EXT_VGLUT-GAL4_RigA_20210305T083721'};
%large oval insert 1LO
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T211624_BLUE'};
% ragged insert 
% expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/singlecolormarkers_try2/pilot24RGB_JHS_K_85321_GtACR1_RigA_20210313T212914_BLUE-oval'};

expis = 1;
% platebowls = {};

% for expi = 1:numel(expdirs),
%   expdir = expdirs{expi};
%   metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%   platebowl = sprintf('%d%s',metadata.plate,metadata.bowl);
%   if ~ismember(platebowl,platebowls),
%     expis(end+1) = expi; %#ok<SAGROW>
%     platebowls{end+1} = platebowl; %#ok<SAGROW>
%   end
% end
% 
% expi = expis(1);
% for ming (no plate)
% expis = expi;
%%

% for i = 1:numel(expis),
%   
%   expi = expis(i);
%   
%   expdir = expdirs{expi};
%   
%   moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
%   annfilename = fullfile(expdir,dataloc_params.annfilestr);
%   [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
%   mu = read_ann(annfilename,'background_median');
%   metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%   
%   nr = headerinfo.nr;
%   nc = headerinfo.nc;
%   mu = reshape(mu,[nc,nr]);
%   
%   figure(i);
%   clf;
%   imagesc(mu); axis image;
%   
% end

%%
for i = 1:numel(expis),
    expi = expis(i);
    expdir = expdirs{expi};
    
    moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
    annfilename = fullfile(expdir,dataloc_params.annfilestr);
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


imwrite(uint8(template),'RegistraionMark_FlyBowlLEDON.png');

%% distance from the corner
x = mean(xlim);
y = mean(ylim);
nc = headerinfo.nc;
nr = headerinfo.nr;
corners = [1,1,nc,nc;1,nr,nr,1];
dcorner = sqrt(min((corners(1,:)-x).^2 + (corners(2,:)-y).^2));
imr = min(nr,nc);

dcornerfrac = dcorner / imr

% maxDistCornerFrac_BowlLabel = .12;