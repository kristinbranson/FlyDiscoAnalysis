function PlotTimeSeries(obj,fly,fnsplot,varargin)

[hfig,hax,tstart,tend] = myparse(varargin,'hfig',[],'hax',[],'tstart',1,'tend',inf);
axparams = {numel(fnsplot),1,[.07,.01;.05,.01]};
[hax,hfig] = get_axes(hax,hfig,'axparams',axparams);

tstart = max(1,tstart);
tend = min(tend,obj.trx(fly).nframes);
t0 = min([obj.trx(obj.movie2flies{obj.fly2movie(fly)}).timestamps]);
for i = 1:numel(fnsplot),
  fn = fnsplot{i};
  n = size(obj.trx(fly).(fn),2);
  t = obj.trx(fly).timestamps;
  t = t - t0;
  if obj.trx(fly).nframes > n,
    d = obj.trx(fly).nframes - n;
    t = conv(t,ones(1,d+1),'valid');
  end
  idx = t >= tstart & t <= tend;
  plot(hax(i),t(idx),obj.trx(fly).(fn)(idx),'k.-');
  ylabel(hax(i),fn,'interpreter','none');
  axisalmosttight([],hax(i));
  if i ~= numel(fnsplot),
    set(hax(i),'XTickLabel',{});
  end
end
xlabel(hax(end),'Time (s)');
set(hax,'XLim',[tstart-1,tend+1],'ticklength',[.005,.0051]);
linkaxes(hax,'x');