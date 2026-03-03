function [him] = PlotAnatomyImageOnAxes(line_name,x,y,w,hax,hdata,varargin)

him = [];
anatomydir = '/nobackup/branson/AverageAnatomyData20130618';

[anatomydir] = myparse(varargin,...
  'anatomydir',anatomydir);

filename = fullfile(anatomydir,sprintf('meanim_%s.png',line_name));
if exist(filename,'file'),
  im = imread(filename);
else
  warning('No image available for line %s',line_name);
  return;
end

isholdval = ishold(hax);
hold(hax,'on');
[nr,nc,~] = size(im);
h = w*nr/nc;

him = image(x+[-w/2,w/2],y+[-h/2,h/2],im,'Parent',hax);
hchil = get(hax,'Children');
hchil(hchil==him) = [];
[ism,idx] = ismember(hdata,hchil);
if ~all(ism),
  warning('Input data handle not a child of the axes');
else
  
  i = max(idx);
  hchil = [hchil(1:i);him;hchil(i+1:end)];
  set(hax,'Children',hchil);
  
end

%axis(hax,'image');

if ~isholdval,
  hold(hax,'off');
end