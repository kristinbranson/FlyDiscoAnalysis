%% set up path

rootdatadir = '/groups/branson/bransonlab/projects/olympiad/WingTrackedData';
addpath ../misc;
addpath ../filehandling;
addpath /groups/branson/home/bransonk/olympiad_bransonlab/FlyBowlAnalysis;

%% parameters

moviefile = '/groups/branson/bransonlab/example_movies/AVT_camera/20120719_courtship/20120719_shelby.ufmf';
annfile = '/groups/branson/bransonlab/example_movies/AVT_camera/20120719_courtship/20120719_shelby_allflies.ufmf.ann';
bg_nframes = 200;


%% read stuff in

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);


%% compute background model

if ~isempty(annfile),
  bgmed = read_ann(annfile,'background_center');
  bgmed = reshape(bgmed,[headerinfo.nc,headerinfo.nr])';
else
  buffer = readframe(1);
  buffer = repmat(buffer,[1,1,bg_nframes]);
  frames = round(linspace(1,nframes,bg_nframes));
  for i = 1:bg_nframes,
    t = frames(i);
    buffer(:,:,i) = readframe(t);
  end
  
  bgmed = double(median(single(buffer),3));
end

%% get regions of interest

hfig = 1;
figure(hfig);
clf;
imagesc(readframe(round(1+nframes)/2));
axis image;
colormap gray;
title('Click to add next region of interest.');

while true,
  

  
  
end