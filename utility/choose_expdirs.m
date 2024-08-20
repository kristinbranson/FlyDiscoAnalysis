function [experiments_chosen,ngal4,ncontrol,idx] = choose_expdirs(experiments,ngal4,ncontrol,varargin)

[doseparatecontrols,weight_order] = myparse(varargin,'doseparatecontrols',true,...
  'weight_order',{'date','genotype','bowl','rig'});

control_line = 'pBDPGAL4U';

if doseparatecontrols,

  % split into control and gal4 lines
  iscontrol = strcmp(control_line,{experiments.line});
else
  iscontrol = false(1,numel(experiments));
end
ngal4 = min(nnz(~iscontrol),ngal4);
ncontrol = min(nnz(iscontrol),ncontrol);
idxcontrol = find(iscontrol);
idxgal4 = find(~iscontrol);

idx = choose_expdirs_helper(experiments(iscontrol),ncontrol,weight_order);
idxcontrol = idxcontrol(idx);
idx = choose_expdirs_helper(experiments(~iscontrol),ngal4,weight_order);
idxgal4 = idxgal4(idx);
idx = [idxgal4,idxcontrol];
experiments_chosen = experiments(idx);

function idxchosen = choose_expdirs_helper(experiments,nout,weight_order)

idxchosen = [];

if nout <= 0,
  return;
end

format = 'yyyymmddTHHMMSS';
dates = datenum({experiments.exp_datetime},format)';

nexpdirs = length(experiments);
rigs = [experiments.rig];
bowls = {experiments.bowl};
genotypes = cell(1,numel(experiments));
for i = 1:numel(experiments),
  genotypes{i} = sprintf('%s__%s',experiments(i).line,experiments(i).effector);
end
maxdate = max(dates) - min(dates) + 1;

isselected = false(1,nexpdirs);
nselected = 0;
nmatches_rig = zeros(1,nexpdirs);
nmatches_bowl = zeros(1,nexpdirs);
nmatches_genotype = zeros(1,nexpdirs);
mind_date = inf(1,nexpdirs);

% choose the first expdir as the last-collected experiment
[~,expdiri] = max(dates);

while true,
  
  % output selected directory
  fprintf('Selecting expdiri = %d, rig = %d, bowl = %s, genotype = %s\n',expdiri,rigs(expdiri),bowls{expdiri},genotypes{expdiri});
  idxchosen(end+1) = expdiri; %#ok<AGROW>
  
  % set that this directory is selected
  isselected(expdiri) = true;
  nselected = nselected + 1;
  
  % break if we've selected enough
  if nselected >= nout || all(isselected),
    break;
  end

  % check for bowl, genotype matches
  idx = rigs(expdiri)==rigs;
  nmatches_rig(idx) = nmatches_rig(idx) + 1;
  idx = strcmpi(bowls{expdiri},bowls);
  nmatches_bowl(idx) = nmatches_bowl(idx) + 1;
  idx = strcmpi(genotypes{expdiri},genotypes);
  nmatches_genotype(idx) = nmatches_genotype(idx) + 1;
  
  % compute distance in date
  mind_date = min(mind_date,abs(dates(expdiri)-dates));

  % combine
  factor = 1;
  dscalar = 0;
  for i = 1:numel(weight_order),
    switch weight_order{i},
      case 'date'
        dscalar = dscalar + factor*mind_date;
        factor = factor * maxdate;
      case 'genotype',
        dscalar = dscalar + factor*(nexpdirs-nmatches_genotype);
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
[genotypesselected,~,j] = unique(genotypes(isselected));
nselectedpergenotype = hist(j,1:length(genotypesselected));

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
fprintf('nselected per genotype:\n');
for i = 1:length(genotypesselected),
  fprintf('  [%s: %d]\n',genotypesselected{i},nselectedpergenotype(i));
end
fprintf('\n');