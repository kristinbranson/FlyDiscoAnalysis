%% set up path

%rootdatadir = '/groups/branson/bransonlab/projects/olympiad/WingTrackedData';
addpath ../misc;
addpath ../filehandling;
addpath ../../Ctrax/matlab/netlab;
%addpath /groups/branson/home/bransonk/olympiad_bransonlab/FlyBowlAnalysis;
rootdatadir = '/groups/branson/bransonlab/projects/FlyBowl/Courtship_plates_movies';
experiment_names = {'CantonS_TrpA_Rig2Plate17BowlC_20121210T150843'
  'CantonS_TrpA_Rig2Plate17BowlD_20121210T150854'
  'CantonS_TrpA_Rig2Plate17BowlA_20121210T150418'
  'CantonS_TrpA_Rig2Plate17BowlB_20121210T150427'};

%% parameters


roifilestr = 'roiinfo.mat';
maxroidcenter = 25;
roiradiusfactor = .975;

roiparams = struct;
roiparams.cannythresh = [.03,.06];
roiparams.cannysigma = 3;
roiparams.maxdcenter = 10;
roiparams.maxdradius = 10;
roiparams.nbinscenter = 21;
roiparams.nbinsradius = 21;
roiparams.meanroiradius = 121.4683;
roiparams.roimus = [
  262.5975  253.3691
  514.7065  256.8469
  766.6785  261.0128
  385.1247  508.2398
  638.3260  511.2702
  889.5829  514.7079
  255.6345  758.9520
  507.6222  762.8746
  760.0276  766.5262
  133.7526  507.6416
  ];

bg_nframes = 20;
estimateroilocs = false;

%% estimate roi locations

if estimateroilocs,

%% read in all roi data

roidata = [];

for expi = 1:numel(experiment_names),

  experiment_name = experiment_names{expi};
  expdir = fullfile(rootdatadir,experiment_name);
  roifilename = fullfile(rootdatadir,experiment_names{expi},roifilestr);
  
  if exist(roifilename,'file'),
    roidata = structappend(roidata,load(roifilename));
    roidata(end).experiment_name = experiment_name;
  end

end

%% get a background median

moviefile = fullfile(rootdatadir,experiment_name,'movie.ufmf');
annfile = fullfile(rootdatadir,experiment_name,'movie.ufmf.ann');

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
imheight = headerinfo.nr;
imwidth = headerinfo.nc;

%% compute background model

buffer = readframe(1);
buffer = repmat(buffer,[1,1,bg_nframes]);
frames = round(linspace(1,nframes,bg_nframes));
for i = 1:bg_nframes,
  t = frames(i);
  buffer(:,:,i) = readframe(t);
end

bgmed = uint8(median(single(buffer),3));

%% plot all rois 

clf;
imagesc(bgmed);
axis image;
colormap gray;
hold on;
colors = lines(numel(roidata));
for i = 1:numel(roidata),
  for j = 1:numel(roidata(i).rois),
    plot(roidata(i).rois{j}(:,1),roidata(i).rois{j}(:,2),'.-','color',colors(i,:));
  end
end

%% greedily match

% fit circles
for i = 1:numel(roidata),
  roidata(i).centerx = nan(1,numel(roidata(i).rois));
  roidata(i).centery = nan(1,numel(roidata(i).rois));
  roidata(i).radius = nan(1,numel(roidata(i).rois));
  for j = 1:numel(roidata(i).rois),
    [roidata(i).centerx(j),roidata(i).centery(j),roidata(i).radius(j)] = ...
      fit_circle_to_points(roidata(i).rois{j}(:,1),roidata(i).rois{j}(:,2));
  end
end

clf;
imagesc(bgmed);
axis image;
colormap gray;
hold on;
colors = lines(numel(roidata));
thetas = linspace(0,2*pi,100);
for i = 1:numel(roidata),
  for j = 1:numel(roidata(i).rois),
    plot(roidata(i).centerx(j)+roiradiusfactor*roidata(i).radius(j)*cos(thetas),...
      roidata(i).centery(j)+roiradiusfactor*roidata(i).radius(j)*sin(thetas),...
      '.-','color',colors(i,:));
  end
end


roimus = zeros(0,2);
roiradii = [];
roiidx = false(numel(roidata),0);
roins = [];
for i = 1:numel(roidata),
  for j = 1:numel(roidata(i).rois),
    idxcurr = find(~roiidx(i,:));
    if ~isempty(idxcurr),
      [d,kk] = min(sqrt( (roimus(idxcurr,1)-roidata(i).centerx(j)).^2 + (roimus(idxcurr,2)-roidata(i).centery(j)).^2 ));
      k = idxcurr(kk);
    else
      d = inf;
    end
    if d > maxroidcenter,
      k = size(roimus,1)+1;
      roins(k) = 0;
      roimus(k,:) = 0;
      roiradii(k) = 0;
    end
    
    roiidx(i,k) = true;
    roimus(k,:) = (roimus(k,:)*roins(k) + [roidata(i).centerx(j),roidata(i).centery(j)]) / (roins(k)+1);
    roiradii(k) = (roiradii(k)*roins(k) + roidata(i).radius(j)) / (roins(k)+1);
    roins(k) = roins(k) + 1;
  end
end

plot(roimus(:,1),roimus(:,2),'kx');

meanroiradius = nanmean(roiradii)*roiradiusfactor;

%% create parameters

roiparams.meanroiradius = meanroiradius;
roiparams.roimus = roimus;

end

%% detect

all_bgmed = cell(1,numel(experiment_names));
all_roicenterx = cell(1,numel(experiment_names));
all_roicentery = cell(1,numel(experiment_names));
all_roiradii = cell(1,numel(experiment_names));
all_roiscores = cell(1,numel(experiment_names));

for expi = 1:numel(experiment_names),
  experiment_name = experiment_names{expi};

  moviefile = fullfile(rootdatadir,experiment_name,'movie.ufmf');

  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
  imheight = headerinfo.nr;
  imwidth = headerinfo.nc;

  % compute background model
  fprintf('Computing background model for %s...\n',experiment_name);
  buffer = readframe(1);
  buffer = repmat(buffer,[1,1,bg_nframes]);
  frames = round(linspace(1,nframes,bg_nframes));
  for i = 1:bg_nframes,
    t = frames(i);
    buffer(:,:,i) = readframe(t);
  end
  
  bgmed = uint8(median(single(buffer),3));
  
  % detect
  
  fprintf('Detecting rois...\n');

  [roicenterx,roicentery,roiradii,roiscores] = detectCourtshipBowlROIs(bgmed,roiparams);
  
  clf;
  imagesc(bgmed);
  hold on;
  axis image;
  colormap gray;
  thetas = linspace(0,2*pi,100);
  for i = 1:size(roimus,1),
    plot(roicenterx(i)+roiradii(i)*cos(thetas),...
      roicentery(i)+roiradii(i)*sin(thetas),'r-');
  end
  
  drawnow;
  
  % store
  
  all_bgmed{expi} = bgmed;
  all_roicenterx{expi} = roicenterx;
  all_roicentery{expi} = roicentery;
  all_roiradii{expi} = roiradii;
  all_roiscores{expi} = roiscores;
  
  
end

