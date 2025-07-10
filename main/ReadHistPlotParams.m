function hist_plot_params = ReadHistPlotParams(histplotparamsfile)

specialfns = {'nbinscombine'};
plotstyles = {'log','linear'};

params = ReadParams(histplotparamsfile);

hist_plot_params = struct;
hist_plot_params.groups = [];

groupfns = setdiff(fieldnames(params),specialfns);
for i = 1:numel(groupfns),
  groupfn = groupfns{i};
  plotstyle = params.(groupfn){1};
  assert(ismember(plotstyle,plotstyles));
  group = struct('name',groupfn,'plotstyle',plotstyle,'features',{params.(groupfn)(2:end)});
  hist_plot_params.groups = structappend(hist_plot_params.groups,group);
end

specialfns1 = intersect(specialfns,fieldnames(params));
for i = 1:numel(specialfns1),
  hist_plot_params.(specialfns1{i}) = params.(specialfns1{i});
end