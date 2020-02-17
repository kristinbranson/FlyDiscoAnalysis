 
close all;
figure(1);
clf;
%set(1,'Position',[252   488   598   615]);
axes('Position',[0,0,1,1]);
im = readframe(1);
him = imagesc(im,[0,255]);
colormap gray;
htext = text(0,0,'00:00:00','Color','w','HorizontalAlignment','Left','VerticalAlignment','top','FontSize',30);
truesize;
% ax = [95.3181  535.3794   11.4318  464.0032];
%ax = [104.75 1654.25 109.218937446444 1664.27463581834];
axis off;
%axis(ax);
axis equal;

%axis(ax);
%set(htext,'Position',[ax(1),ax(3),0],'Color','w')

%%
for interval = 1:numel(firstframes),
for i = 1:maxnframes,
t = firstframes(interval)+i-1;
im = readframe(t);
set(him,'CData',im);
set(htext,'String',datestr((headerinfo.timestamps(t)-headerinfo.timestamps(firstframes(1)))/24/3600,'HH:MM:SS'));
fr = getframe(1);
writeVideo(aviobj,fr);
end
end