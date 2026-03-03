function [outmeanfile,outpvaluefile,outmatfile,meanavifile,pvalueavifile] = SaveAverageAnatomyImage(meanim,pvalue,outbasename,varargin)

[clim_pvalue,clim_meanim,cm_meanim,cm_pvalue] = myparse(varargin,'clim_pvalue',[0,.01],...
  'clim_meanim',[0,nan],'cm_meanim',kjetsmooth(256),'cm_pvalue',flipud(kjetsmooth(256)));

outmeanfile = [outbasename,'_meanim.png'];
Irgb = im2uint8(colormap_image(max(meanim,[],3)',cm_meanim,clim_meanim));
imwrite(Irgb,outmeanfile,'png');

outpvaluefile = [outbasename,'_pvalue.png'];
Irgb = im2uint8(colormap_image(min(pvalue,[],3)',cm_pvalue,clim_pvalue));
imwrite(Irgb,outpvaluefile,'png');

outmatfile = [outbasename,'.mat'];
save('-v7.3',outmatfile,'meanim','pvalue','clim_pvalue','clim_meanim','cm_pvalue','cm_meanim');

%% make a z-stack movie

hfig = gcf;
clf;
set(hfig,'Units','pixels','Position',[10,10,1024,512]);
hax = axes('Position',[0,0,1,1]);
z = 1;
him = imagesc(meanim(:,:,1)');
colormap(cm_meanim);
if isnan(clim_meanim(1)),
  clim_mean(1) = 0;
end
if isnan(clim_meanim(2)),
  clim_mean(2) = max(meanim(:));
end
set(gca,'CLim',clim_mean);
hcb = colorbar('Location','East');
set(hcb,'XColor','w','YColor','w');
axis image off;
htext = text(0,0,sprintf(' z = %d',z),'Color','w','HorizontalAlignment','left',...
  'VerticalAlignment','top','FontSize',24);
set(hfig,'Visible','off');

meanuavifile = [outbasename,'_meanstack.avi'];
meanavifile = [outbasename,'_meanstack_xvid.avi'];
aviobj = VideoWriter(meanuavifile,'Uncompressed AVI');
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
  meanuavifile,meanavifile,newwidth,newheight);
status = system(cmd);
if status ~= 0,
  warning('Could not compress video %s',meanuavifile);
else
  delete(meanuavifile);
end

%% make a z-stack movie for the p-value

hfig = gcf;
clf;
set(hfig,'Units','pixels','Position',[10,10,1024,512]);
hax = axes('Position',[0,0,1,1]);
z = 1;
him = imagesc(pvalue(:,:,1)');
colormap(cm_pvalue);
if isnan(clim_pvalue(1)),
  clim_pvalue(1) = 0;
end
if isnan(clim_pvalue(2)),
  clim_pvalue(2) = max(pvalue(:));
end
set(gca,'CLim',clim_pvalue);
hcb = colorbar('Location','East');
poscb = get(hcb,'Position');
poscb(end) = poscb(end)*.95;
poscb(1) = poscb(1)*.99;
set(hcb,'Position',poscb);
set(hcb,'XColor','w','YColor','w');
axis image off;
htext = text(0,0,sprintf(' z = %d',z),'Color','w','HorizontalAlignment','left',...
  'VerticalAlignment','top','FontSize',24);
set(hfig,'Visible','off');

pvalueuavifile = [outbasename,'_pvaluestack.avi'];
pvalueavifile = [outbasename,'_pvaluestack_xvid.avi'];
aviobj = VideoWriter(pvalueuavifile,'Uncompressed AVI');
aviobj.FrameRate = 5;
open(aviobj);

fr = getframe_invisible(hax);
[height,width,~] = size(fr);
fprintf('Size of frame is %d x %d\n',height,width);
gfdata = getframe_initialize(hax);
[fr,height,width] = getframe_invisible_nocheck(gfdata,[height,width],false,false);
writeVideo(aviobj,fr);

for z = 2:size(pvalue,3),
  set(him,'CData',pvalue(:,:,z)');
  set(htext,'String',sprintf(' z = %d',z));
  drawnow;
  fr = getframe_invisible_nocheck(gfdata,[height,width],false);
  writeVideo(aviobj,fr);
end
close(aviobj);

newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d -msglevel all=2',...
  pvalueuavifile,pvalueavifile,newwidth,newheight);
status = system(cmd);
if status ~= 0,
  warning('Could not compress video %s',pvalueuavifile);
else
  delete(pvalueuavifile);
end
