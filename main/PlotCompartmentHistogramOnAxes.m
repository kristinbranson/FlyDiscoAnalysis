function [hpatch] = PlotCompartmentHistogramOnAxes(line_name,x,y,wfrac,hfrac,hax,hdata,int_manual,scalecolorby,varargin)

hpatch = [];
persistent colors;
if isempty(colors),
  colors = jet(256);
end

compartments = setdiff(fieldnames(int_manual),'line_names');
[compartments] = myparse(varargin,'orderedcompartments',compartments);

idx = find(strcmp(line_name,int_manual.line_names),1);
if isempty(idx),
  warning('No annotation available for line %s',line_name);
  return;
end

isholdval = ishold(hax);
hold(hax,'on');

ncompartments = numel(compartments);

datacurr = nan(1,ncompartments);
for i = 1:ncompartments,
  datacurr(i) = int_manual.(compartments{i})(idx);
end
datacurr(isnan(datacurr)) = 0;

%colori = randsample(size(colors,1),1);

colori = min(size(colors,1),max(1,round(nansum(datacurr)/scalecolorby*size(colors,1))));
color = colors(colori,:);

ax = axis(hax);
w = diff(ax(1:2))*wfrac;
h = diff(ax(3:4))*hfrac;

hpatch = nan(1,2);
hpatch(1) = plot(x+w/2*[-1,-1,1],y+h/2*[1,-1,-1],'k-');
htext = [];

xplot = (x-w/2)+bsxfun(@plus,1:ncompartments,[-.5;.5])/ncompartments*w;
yplot = (y-h/2)+[zeros(1,ncompartments);datacurr]/5*h;
xplot = xplot([1,1,2,2],:);
xplot = [xplot(:);xplot(1)];
yplot = yplot([1,2,2,1],:);
yplot = [yplot(:);yplot(1)];
xplot(4:4:end-5) = [];
yplot(4:4:end-5) = [];
hpatch(2) = patch(xplot,yplot,color,'EdgeColor','none','Parent',hax);

for i = find(datacurr > 2),
  xcurr = (x-w/2)+i/ncompartments*w;
  ycurr = (y-h/2)+datacurr(i)/5*h;
  htext(end+1) = text(xcurr,ycurr,compartments{i},'HorizontalAlignment','center','VerticalAlignment','bottom','Parent',hax); %#ok<AGROW>
end
% i = maxi;
% xcurr = (x-w/2)+i/ncompartments*w;
% ycurr = (y-h/2)+int_manual.(compartments{i})(idx)/5*h;
% hpatch(ncompartments+2) = text(xcurr,ycurr,compartments{i},'HorizontalAlignment','center','VerticalAlignment','bottom');
hpatch = [hpatch,htext];

hchil = get(hax,'Children');
hchil(ismember(hchil,hpatch)) = [];
[ism,idx] = ismember(hdata,hchil);
if ~all(ism),
  warning('Input data handle not a child of the axes');
else
  
  i = max(idx);
  hchil = [hchil(1:i);hpatch';hchil(i+1:end)];
  set(hax,'Children',hchil);
  
end

if ~isholdval,
  hold(hax,'off');
end