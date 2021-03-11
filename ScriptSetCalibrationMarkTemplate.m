%% set up paths
    
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
% settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings/';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings/';


%% parameters

% % analysis_protocol = '20150717_flybubble_flybowlMing';
analysis_protocol = '20210219_flybubble_dickson_RGBGtACR1testing';
% expfile = '/groups/branson/home/robiea/Projects_data/Ming/explist_forregistration';
expfile = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/FlyBubbleRGB/locomotionGtACR1_24_RGB_EXT_VGLUT-GAL4_RigA_20210305T083721'
datalocparamsfilestr = 'dataloc_params.txt';
dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));

expdirs = importdata(expfile);
expdirs(cellfun(@isempty,expdirs)) = [];

expis = [];
platebowls = {};

for expi = 1:numel(expdirs),
  expdir = expdirs{expi};
  metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
  platebowl = sprintf('%d%s',metadata.plate,metadata.bowl);
  if ~ismember(platebowl,platebowls),
    expis(end+1) = expi; %#ok<SAGROW>
    platebowls{end+1} = platebowl; %#ok<SAGROW>
  end
end

expi = expis(1);

%%

for i = 1:numel(expis),
  
  expi = expis(i);
  
  expdir = expdirs{expi};
  
  moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
  annfilename = fullfile(expdir,dataloc_params.annfilestr);
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefilename);
  mu = read_ann(annfilename,'background_median');
  metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
  
  nr = headerinfo.nr;
  nc = headerinfo.nc;
  mu = reshape(mu,[nc,nr]);
  
  figure(i);
  clf;
  imagesc(mu); axis image;
  
end

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


imwrite(uint8(template),'CirleTemplate_FlyBowlMing.png');

%% distance from the corner
x = mean(xlim);
y = mean(ylim);
corners = [1,1,nc,nc;1,nr,nr,1];
dcorner = sqrt(min((corners(1,:)-x).^2 + (corners(2,:)-y).^2));
imr = min(nr,nc);
dcornerfrac = dcorner / imr;

maxDistCornerFrac_BowlLabel = .12;