function [fracsmaller_stat,fracbigger_stat] = ComputePValueBySampling(stati,filename,outfiledir,iscircstat)

rng default;

if ischar(stati),
  stati = str2double(stati);
end

load(filename);

% Non-control lines to evaluate. Position-agnostic: control lines may be
% anywhere and there may be more than one (e.g. VNC has YNA + JHS). Output is
% indexed compactly (1:nlines_noncontrol) and the caller maps it back via
% noncontrol_line_idx. Fall back to "all but last" for old outmatfiles.
if ~exist('noncontrol_line_idx','var') || isempty(noncontrol_line_idx)
  noncontrol_line_idx = 1:nlines-1;
end
if ~exist('nlines_noncontrol','var') || isempty(nlines_noncontrol)
  nlines_noncontrol = numel(noncontrol_line_idx);
end

fracsmaller_stat = nan(nlines_noncontrol,1);
fracbigger_stat = nan(nlines_noncontrol,1);

statfn = statfns{stati}; %#ok<USENS>
% statfn is the full (display) name; statkey is the MATLAB-safe struct field key
% (== statfn unless the full name exceeds 63 char). Fall back to statfn for
% older outmatfiles saved without statkeys.
if exist('statkeys','var') && ~isempty(statkeys)
  statkey = statkeys{stati};
else
  statkey = statfn;
end

fprintf('Statistic %s: %d / %d\n',statfn,stati,nstats);

% --- slow-control fallback (read from outmatfile; back-compatible defaults) ---
% isnanstat: too sparse for any null -> return all NaN.
% usealldata: build the null from ALL good sets (population) instead of control.
% mingoodexps: a set is a valid resampling unit only with this many good exps.
if exist('isnanstat','var') && ~isempty(isnanstat) && isnanstat(stati),
  fprintf('Stat %s: too sparse for a null (isnanstat), returning NaN\n',statfn);
  if ~exist(outfiledir,'dir'), mkdir(outfiledir); end
  save(fullfile(outfiledir,sprintf('PvaluesForStat%03d.mat',stati)),'fracsmaller_stat','fracbigger_stat');
  return;
end
use_all = exist('usealldata','var') && ~isempty(usealldata) && usealldata(stati);
if exist('mingoodexpsperset','var') && ~isempty(mingoodexpsperset),
  mingoodexps = mingoodexpsperset;
else
  mingoodexps = 2;
end

isgoodset = ~isinf(setstats.normmeans.(statkey)) & ...
  ~isnan(setstats.normmeans.(statkey)) & ...
  setstats.nexps.(statkey) >= minnexps;
isgoodexp = ~isinf(allstats.(statkey)) & ...
  ~isnan(allstats.(statkey)) & ...
  allnframestotal.(statkey) >= minnframes(stati);

% resampling pool: control sets, or ALL good sets under the fallback
if use_all,
  isgoodpoolset = isgoodset;
else
  isgoodpoolset = isgoodset & setiscontrol;
end
if ~any(isgoodpoolset),
  fprintf('No good sets in pool for statistic %s, not computing p-values\n',statfn);
  return;
end


idxgoodpoolset = find(isgoodpoolset);

set2expcurr = nan(nsets,maxnexps);
ngoodpoolset = nan(1,nsets);
for seti = idxgoodpoolset,
  tmp = find(exp2setidx==seti & isgoodexp);
  ngoodpoolset(seti) = numel(tmp); % added AR 20260113 fix bug of NaN in setsampleis
  set2expcurr(seti,1:numel(tmp)) = tmp;
end 

isgoodline = ~isnan(linestats.normmeans.(statkey)) & ...
  ~isinf(linestats.normmeans.(statkey));

for ii = 1:nlines_noncontrol,
  linei = noncontrol_line_idx(ii);   % actual line index (controls excluded, any position)

  if mod(ii,100) == 1,
    fprintf('Stat %s (%d / %d), line %s (%d / %d)\n',statfn,stati,nstats,...
      linestats.line_names{linei},ii,nlines_noncontrol);
  end

  if ~isgoodline(linei),
    fracsmaller_stat(ii) = 1;
    fracbigger_stat(ii) = 1;
    continue;
  end
  
  % number of experiments in each set
  setidxcurr = find(set2lineidx==linei & isgoodset);
  if isempty(setidxcurr),
    %fprintf('No good sets found for line %s, skipping.\n',linestats.line_names{linei});
    continue;
  end
  
  nexpscurr = setstats.nexps.(statkey)(setidxcurr);
  nexpscurr = sort(nexpscurr,'descend');
  nsetscurr = numel(nexpscurr);
  
  % sample
  setnormmu = nan(nsamples,nsetscurr);
  skipline = false;
  for setii = 1:nsetscurr,
    setidxallowed = find(isgoodpoolset & ...
      setstats.nexps.(statkey) >= nexpscurr(setii) & ...
      ngoodpoolset >= nexpscurr(setii)); % added AR 20260113 fix bug of NaN in setsampleis
    nsub = 1; % defensive: keep defined if the loop below does not run (nexpscurr<=2)
    for nsub = 1:nexpscurr(setii)-mingoodexps,   % floor nexpscurr at mingoodexps (=2): never resample <2 exps/set
      if numel(setidxallowed) > 1,
        break;
      end
      % keep the ngoodpoolset constraint in the fallback too, else
      % set2expcurr can contain NaN for sparse (binned) stats -> bad index
      setidxallowed = find(isgoodpoolset & ...
        setstats.nexps.(statkey) >= nexpscurr(setii)-nsub & ...
        ngoodpoolset >= nexpscurr(setii)-nsub);
    end
    nsub = nsub-1;
    if nsub > 0,
      fprintf('For stat %s, line %s, set %d, needed to consider sets %d exps smaller than than this set\n',statfn,linestats.line_names{linei},setii,nsub);
      nexpscurr(setii) = nexpscurr(setii)-nsub;
    end
    % Re-select control sets consistent with the FINAL nexpscurr(setii). The
    % fallback loop can settle on a nexpscurr whose threshold does not match the
    % last setidxallowed it kept (off-by-one on the exhaustion path), letting a
    % set with too few good exps slip in -> set2expcurr NaN -> crash. Re-finding
    % here guarantees every chosen set has >= nexpscurr(setii) good experiments.
    setidxallowed = find(isgoodpoolset & ...
      setstats.nexps.(statkey) >= nexpscurr(setii) & ...
      ngoodpoolset >= nexpscurr(setii));
    if isempty(setidxallowed),
      % no control set has enough frame-good experiments for this line/set:
      % this line cannot be sampled for this (sparse) stat -> leave p-value NaN
      fprintf('For stat %s, line %s: no control set with enough good exps; line p-value left NaN\n',statfn,linestats.line_names{linei});
      skipline = true;
      break;
    end
    setsampleis = randsample(setidxallowed,nsamples,true);
    
    % choose the experiments per set
    % controlnexpscurr = setstats.nexps.(statkey)(setsampleis); % added AR 20260113 fix bug of NaN in setsampleis
    controlnexpscurr = ngoodpoolset(setsampleis); % added AR 20260113 fix bug of NaN in setsampleis
    % by default, use the first experiments
    expsampleis = set2expcurr(setsampleis,1:nexpscurr(setii));
    
    % for sets with more experiments, sample without replacement
    tmpidx = find(controlnexpscurr > nexpscurr(setii));
    for tmpi = tmpidx,
      expsampleis(tmpi,:) = set2expcurr(setsampleis(tmpi),randsample(controlnexpscurr(tmpi),nexpscurr(setii)));
    end
    
    tmp = allstats.(statkey)(expsampleis);
    tmp = reshape(tmp,[nsamples,nexpscurr(setii)]);
    if iscircstat(stati)
        mu = circ_mean(tmp,[],2);
        % setnormmu(:,setii) = mu - setstats.controlmeans.(statkey)(setsampleis)' + controlmean(stati);
        % Circular offset
        offset = circ_dist(mu, setstats.controlmeans.(statkey)(setsampleis)');
        % Add offset to control mean (circular addition = regular addition + wrapping)
        setnormmu(:, setii) = angle(exp(1i * (controlmean(stati) + offset)));    
    else
    mu = mean(tmp,2);
    setnormmu(:,setii) = mu - setstats.controlmeans.(statkey)(setsampleis)' + controlmean(stati);
    end
    
  end
  if skipline,
    continue;   % unsamplable line for this sparse stat -> p-value stays NaN
  end
  if iscircstat(stati)
      mu = circ_mean(setnormmu,[],2);
      % circ_dist(line,mu) > 0 means line > mu, so fraction of controls smaller than line
      fracsmaller_stat(ii) = nnz(0<circ_dist(linestats.normmeans.(statkey)(linei),mu))/nsamples;
      fracbigger_stat(ii) = nnz(0>circ_dist(linestats.normmeans.(statkey)(linei),mu))/nsamples;
  else
      mu = mean(setnormmu,2);
      fracsmaller_stat(ii) = nnz(mu<linestats.normmeans.(statkey)(linei))/nsamples;
      fracbigger_stat(ii) = nnz(mu>linestats.normmeans.(statkey)(linei))/nsamples;
  end
end

if ~exist(outfiledir,'dir'),
  mkdir(outfiledir);
end

outfilename = fullfile(outfiledir,sprintf('PvaluesForStat%03d.mat',stati));
save(outfilename,'fracsmaller_stat','fracbigger_stat');
