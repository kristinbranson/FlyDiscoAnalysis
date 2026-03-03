function [experiments_chosen,nchosen,idxchosen] = ChooseExperimentsCondition(experiments,nchoose,varargin)

[weight_order] = myparse(varargin,'weight_order',{'screen_type','screen_reason','rig','bowl','date'});

nchoose = min(numel(experiments),nchoose);

idxchosen = [];

if nchoose <= 0,
  return;
end

format = 'yyyymmddTHHMMSS';
dates = datenum({experiments.exp_datetime},format)';

nexpdirs = length(experiments);
screen_types = {experiments.screen_type};
screen_reasons = {experiments.screen_reason};
rigs = [experiments.rig];
bowls = {experiments.bowl};
maxdate = max(dates) - min(dates) + 1;

isselected = false(1,nexpdirs);
nselected = 0;
nmatches_type = zeros(1,nexpdirs);
nmatches_reason = zeros(1,nexpdirs);
nmatches_rig = zeros(1,nexpdirs);
nmatches_bowl = zeros(1,nexpdirs);
mind_date = inf(1,nexpdirs);

% choose the first expdir as the last-collected experiment
[~,expdiri] = max(dates);

while true,
  
  % output selected directory
  fprintf('Selecting expdiri = %d, screen_type = %s, screen_reason = %s, rig = %d, bowl = %s\n',...
    expdiri,screen_types{expdiri},screen_reasons{expdiri},rigs(expdiri),bowls{expdiri});
  idxchosen(end+1) = expdiri; %#ok<AGROW>
  
  % set that this directory is selected
  isselected(expdiri) = true;
  nselected = nselected + 1;
  
  % break if we've selected enough
  if nselected >= nchoose || all(isselected),
    break;
  end

  % check for bowl, genotype matches
  idx = rigs(expdiri)==rigs;
  nmatches_rig(idx) = nmatches_rig(idx) + 1;
  idx = strcmpi(bowls{expdiri},bowls);
  nmatches_bowl(idx) = nmatches_bowl(idx) + 1;
  idx = strcmpi(screen_types{expdiri},screen_types);
  nmatches_type(idx) = nmatches_type(idx) + 1;
  idx = strcmpi(screen_reasons{expdiri},screen_reasons);
  nmatches_reason(idx) = nmatches_reason(idx) + 1;
  
  % compute distance in date
  mind_date = min(mind_date,abs(dates(expdiri)-dates));

  % combine
  factor = 1;
  dscalar = zeros(1,nexpdirs);
  for i = 1:numel(weight_order),
    switch weight_order{i},
      case 'date'
        dscalar = dscalar + factor*mind_date;
        factor = factor * maxdate;
      case 'screen_type',
        dscalar = dscalar + factor*(nexpdirs-nmatches_type);
        factor = factor * nexpdirs;
      case 'screen_reason',
        dscalar = dscalar + factor*(nexpdirs-nmatches_reason);
        factor = factor * nexpdirs;
      case 'bowl',
        dscalar = dscalar + factor*(nexpdirs-nmatches_bowl);
        factor = factor * nexpdirs;
      case 'rig',
        dscalar = dscalar + factor*(nexpdirs-nmatches_rig);
        factor = factor * nexpdirs;
    end        
  end
  %dscalar = mind_date + maxdate*(nexpdirs-nmatches_genotype) + maxdate*nexpdirs*(nexpdirs-nmatches_bowl) + maxdate*nexpdirs^2*(nexpdirs-nmatches_rig);
  dscalar(idxchosen) = -inf;
  
  % choose the farthest
  [~,expdiri] = max(dscalar);
  if isselected(expdiri), 
    error('Sanity check: selecting previously selected expdir'); 
  end
    
end

% statistics
[rigsselected,~,j] = unique(rigs(isselected));
nselectedperrig = hist(j,1:length(rigsselected));
[bowlsselected,~,j] = unique(bowls(isselected));
nselectedperbowl = hist(j,1:length(bowlsselected));
[typesselected,~,j] = unique(screen_types(isselected));
nselectedpertype = hist(j,1:length(typesselected));
[reasonsselected,~,j] = unique(screen_reasons(isselected));
nselectedperreason = hist(j,1:length(reasonsselected));

fprintf('*Summary:\n');
fprintf('nselected = %d\n',nselected);
fprintf('nselected per rig:\n');
for i = 1:length(rigsselected),
  fprintf('  [%d: %d]\n',rigsselected(i),nselectedperrig(i));
end
fprintf('nselected per bowl:\n');
for i = 1:length(bowlsselected),
  fprintf('  [%s: %d]\n',bowlsselected{i},nselectedperbowl(i));
end
fprintf('nselected per screen_type:\n');
for i = 1:length(typesselected),
  fprintf('  [%s: %d]\n',typesselected{i},nselectedpertype(i));
end
fprintf('nselected per screen_reason:\n');
for i = 1:length(reasonsselected),
  fprintf('  [%s: %d]\n',reasonsselected{i},nselectedperreason(i));
end
fprintf('\n');

experiments_chosen = experiments(idxchosen);

nchosen = numel(idxchosen);