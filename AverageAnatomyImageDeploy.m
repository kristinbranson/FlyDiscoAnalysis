function meanim = AverageAnatomyImageDeploy(paramsfile)

% [useallims,maxqi,filsizexy,filsizez,savedir,cm,hfig] = ...
%   myparse(varargin,'useallims',false,'maxqi',.1,...
%   'filsizexy',5,'filsizez',2,'savedir','',...
%   'colormap',kjetsmooth(256),'hfig',[]);

load(paramsfile);

if ~exist(savedir,'dir'),
  mkdir(savedir);
end

datafilenames = {};
for i = 1:numel(line_names),
  filename = fullfile(datadir,sprintf('meanim_%s.mat',line_names{i}));
  if ~exist(filename,'file'),
    fprintf('Line %s is not imaged\n',line_names{i});
    continue;
  end
  datafilenames{end+1} = filename; %#ok<AGROW>
end

nlinescurr = numel(datafilenames);
if nlinescurr == 0,
  error('No lines found\n');
end

nlinescount = 0;
maxvs = [];
minvs = [];
for ii = 1:nlinescurr,

  tmp = load(datafilenames{ii});
  meanimcurr = tmp.meanim / max(tmp.meanim(:));
  
  if nlinescount == 0,
    meanim = meanimcurr;
  else
    meanim = meanim + meanimcurr;
  end

  nlinescount = nlinescount + 1;

  maxvs(nlinescount) = max(meanim(:))/nlinescount;
  minvs(nlinescount) = min(meanim(:))/nlinescount;
  
  outfile = fullfile(savedir,sprintf('meanim%02d.png',nlinescount));
  Irgb = im2uint8(colormap_image(max(meanim/nlinescount,[],3)',cm));
  imwrite(Irgb,outfile,'png');
  
end

meanim = meanim / nlinescount;

outfile = fullfile(savedir,'meanim.png');
Irgb = im2uint8(colormap_image(max(meanim,[],3)',cm));
imwrite(Irgb,outfile,'png');

save(fullfile(savedir,'meanim.mat'),'meanim','datafilenames','maxvs','minvs');

%% make a z-stack movie

hfig = gcf;
clf;
set(hfig,'Units','pixels','Position',[10,10,1024,512]);
hax = axes('Position',[0,0,1,1]);
z = 1;
him = imagesc(meanim(:,:,1)');
colormap(kjetsmooth(256));
set(gca,'CLim',[0,max(meanim(:))]);
hcb = colorbar('Location','East');
set(hcb,'XColor','w','YColor','w');
axis image off;
htext = text(0,0,sprintf(' z = %d',z),'Color','w','HorizontalAlignment','left',...
  'VerticalAlignment','top','FontSize',24);
set(hfig,'Visible','off');

uavifile = fullfile(savedir,'meanstack.avi');
avifile = fullfile(savedir,'meanstack_xvid.avi');
aviobj = VideoWriter(uavifile,'Uncompressed AVI');
aviobj.FrameRate = 5;
open(aviobj);

fr = getframe_invisible(hax);
[height,width,~] = size(fr);
fprintf('Size of frame is %d x %d\n',height,width);
gfdata = getframe_initialize(hax);
[fr,height,width] = getframe_invisible_nocheck(gfdata,[height,width],false,false);
writeVideo(aviobj,fr);

for z = 2:size(meanim,3),
  set(him,'CData',meanim(:,:,z)');
  set(htext,'String',sprintf(' z = %d',z));
  drawnow;
  fr = getframe_invisible_nocheck(gfdata,[height,width],false);
  writeVideo(aviobj,fr);
end
close(aviobj);

newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d -msglevel all=2',...
  uavifile,avifile,newwidth,newheight);
status = system(cmd);
if status ~= 0,
  warning('Could not compress video');
else
  delete(uavifile);
end