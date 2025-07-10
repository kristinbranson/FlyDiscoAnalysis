function [fracsmaller_stat,fracbigger_stat] = ComputePValueBySampling(stati,filename,outfiledir)

rng default;

if ischar(stati),
  stati = str2double(stati);
end

load(filename);

fracsmaller_stat = nan(nlines-1,1);
fracbigger_stat = nan(nlines-1,1);

statfn = statfns{stati}; %#ok<USENS>

fprintf('Statistic %s: %d / %d\n',statfn,stati,nstats);

isgoodset = ~isinf(setstats.normmeans.(statfn)) & ...
  ~isnan(setstats.normmeans.(statfn)) & ...
  setstats.nexps.(statfn) >= minnexps;
isgoodexp = ~isinf(allstats.(statfn)) & ...
  ~isnan(allstats.(statfn)) & ...
  allnframestotal.(statfn) >= minnframes(stati);

isgoodcontrolset = isgoodset & setiscontrol;
if ~any(isgoodcontrolset),
  fprintf('One or less good control sets for statistic %s, not computing p-values\n',statfn);
  return;
end


idxgoodcontrolset = find(isgoodcontrolset);

set2expcurr = nan(nsets,maxnexps);
for seti = idxgoodcontrolset,
  tmp = find(exp2setidx==seti & isgoodexp);
  set2expcurr(seti,1:numel(tmp)) = tmp;
end

isgoodline = ~isnan(linestats.normmeans.(statfn)) & ...
  ~isinf(linestats.normmeans.(statfn));

for linei = 1:nlines-1,
  
  if mod(linei,100) == 1,
    fprintf('Stat %s (%d / %d), line %s (%d / %d)\n',statfn,stati,nstats,...
      linestats.line_names{linei},linei,nlines-1);
  end
  
  %fprintf('Line %s: %d / %d\n',linestats.line_names{linei},linei,nlines);
  
  if ~isgoodline(linei),
    fracsmaller_stat(linei) = 1;
    fracbigger_stat(linei) = 1;
    continue;
  end
  
  % number of experiments in each set
  setidxcurr = find(set2lineidx==linei & isgoodset);
  if isempty(setidxcurr),
    %fprintf('No good sets found for line %s, skipping.\n',linestats.line_names{linei});
    continue;
  end
  
  nexpscurr = setstats.nexps.(statfn)(setidxcurr);
  nexpscurr = sort(nexpscurr,'descend');
  nsetscurr = numel(nexpscurr);
  
  % sample
  setnormmu = nan(nsamples,nsetscurr);
  for setii = 1:nsetscurr,
    
    setidxallowed = find(isgoodcontrolset & ...
      setstats.nexps.(statfn) >= nexpscurr(setii));
    for nsub = 1:nexpscurr(setii)-1,
      if numel(setidxallowed) > 1,
        break;
      end
      setidxallowed = find(isgoodcontrolset & ...
        setstats.nexps.(statfn) >= nexpscurr(setii)-nsub);
    end
    nsub = nsub-1;
    if nsub > 0,
      fprintf('For stat %s, line %s, set %d, needed to consider sets %d exps smaller than than this set\n',statfn,linestats.line_names{linei},setii,nsub);
      nexpscurr(setii) = nexpscurr(setii)-nsub;
    end
    setsampleis = randsample(setidxallowed,nsamples,true);
    
    % choose the experiments per set
    controlnexpscurr = setstats.nexps.(statfn)(setsampleis);
    
    % by default, use the first experiments
    expsampleis = set2expcurr(setsampleis,1:nexpscurr(setii));
    
    % for sets with more experiments, sample without replacement
    tmpidx = find(controlnexpscurr > nexpscurr(setii));
    for tmpi = tmpidx,
      expsampleis(tmpi,:) = set2expcurr(setsampleis(tmpi),randsample(controlnexpscurr(tmpi),nexpscurr(setii)));
    end
    
    tmp = allstats.(statfn)(expsampleis);
    tmp = reshape(tmp,[nsamples,nexpscurr(setii)]);
    mu = mean(tmp,2);
    setnormmu(:,setii) = mu - setstats.controlmeans.(statfn)(setsampleis)';
    
  end
  
  mu = mean(setnormmu,2);
  fracsmaller_stat(linei) = nnz(mu<linestats.normmeans.(statfn)(linei))/nsamples;
  fracbigger_stat(linei) = nnz(mu>linestats.normmeans.(statfn)(linei))/nsamples;
  
end

if ~exist(outfiledir,'dir'),
  mkdir(outfiledir);
end

outfilename = fullfile(outfiledir,sprintf('PvaluesForStat%03d.mat',stati));
save(outfilename,'fracsmaller_stat','fracbigger_stat');
