%% load in data

load CollectedPerFrameStats20130306.mat;
statfns = fieldnames(allstats);
nstats = numel(statfns);

metadata = metadata(idxanalyze);
for i = 1:numel(statfns),
  allstats.(statfns{i})(~idxanalyze) = [];
end
nlines = numel(linestats.line_names);

%% parameters

nsamples = 1000;

%% look at one statistic for control data only

statfn = 'fractime_flyany_framestop';
stati = find(strcmp(statfns,statfn));
linecontrolidx = find(strcmp(linestats.line_names,'pBDPGAL4U'));
expiscontrol = strcmp({metadata.line_name},'pBDPGAL4U');
expidxcontrol = find(expiscontrol);
setiscontrol = set2lineidx == linecontrolidx;
setidxcontrol = find(setiscontrol);
x = allstats.(statfn)(expiscontrol);
setx = setstats.means.(statfn)(setiscontrol);
nexpsperset_control = setnexps(setiscontrol,stati);
[~,exp2set_control] = ismember(setidx(expiscontrol),setidxcontrol);
metadata_control = metadata(expiscontrol);

% weight all experiments the same, and compute the mean
muhat_exp = mean(x);
sigmahat_exp_total = std(x);

% weight all sets the same, and compute the mean
muhat_set = mean(setx(nexpsperset_control>2));
sigmahat_set_total = std(setx(nexpsperset_control>2));

fprintf('Mean %s weighting each experiment the same:    %f\n',statfn,muhat_exp);
fprintf('Mean %s weighting each set w nexps>2 the same: %f\n',statfn,muhat_set);
% These are basically the same, so just use muhat_exp
% Mean fractime_flyany_framestop weighting each experiment the same:    0.417324
% Mean fractime_flyany_framestop weighting each set w nexps>2 the same: 0.417487

fprintf('Std of experiments around the mean:    %f\n',sigmahat_exp_total);
fprintf('Std of sets w nexps>2 around the mean: %f\n',sigmahat_set_total);
% Standard deviation of experiments is bigger than sets
% Std of experiments around the mean:    0.071519
% Std of sets w nexps>2 around the mean: 0.065954

% compute the standard deviation of the experiments around their set means
std_perset = nan(1,numel(setx));
for i = 1:numel(setx),
  idx = find(exp2set_control==i);
  std_perset(i) = sqrt(sum( (x(idx)-setx(i)).^2 )/(numel(idx)-1));
end
% use the mean variance
sigmahat_exp = sqrt(mean(std_perset(nexpsperset_control > 2).^2));

% approximate portion of variance due to differences in sets
mean_nexpsperset = mean(nexpsperset_control(nexpsperset_control > 2));
sigmahat_set = sqrt( sigmahat_set_total.^2 - sigmahat_exp.^2 / mean_nexpsperset );
fprintf('Standard deviation of noise from experiments = %f.\n',sigmahat_exp);
fprintf('Standard deviation of noise from sets = %f.\n',sigmahat_set);
% Standard deviation of noise from experiments = 0.031850.
% Standard deviation of noise from sets = 0.063953.

% how to weight things somewhat properly, at least between exps and sets
n = nexpsperset_control(exp2set_control);
expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
[muhat_weighted,sigma2hat_weighted] = weighted_mean_cov(x(:),expweights(:));
sigmahat_weighted = sqrt(sigma2hat_weighted);

sigmahat_set0 = sigmahat_set;
sigmahat_exp0 = sigmahat_exp;

%% initialize with sigmahat_set and sigmahat_exp, then re-optimize

sigmahat_set_perstat = nan(1,numel(statfns));
sigmahat_exp_perstat = nan(1,numel(statfns));
muhat_weighted_perstat = nan(nlines,numel(statfns));
sigmahat_weighted_perstat = nan(nlines,numel(statfns));

statis = find(~cellfun(@isempty,regexp(statfns,'^fractime.*','once')));

[~,exp2lineidx] = ismember({metadata.line_name},linestats.line_names);

for stati = statis(:)',
  
  statfn = statfns{stati};
  
  x = allstats.(statfn);
  badexpidx = isnan(x);
  setx = setstats.means.(statfn);
  badsetidx = isnan(setx);
  nexpsperset = setnexps(:,stati);
  exp2set = setidx;
  n = nexpsperset(exp2set);
  
  sigmahat_set = sigmahat_set0;
  sigmahat_exp = sigmahat_exp0;
  
  for iter = 1:10,
    
    % given the current values of sigmahat_set and sigmahat_exp, estimate
    % muhat_perline
    expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
    
    sum_perline = accumarray(exp2lineidx(~badexpidx)',x(~badexpidx)'.*expweights(~badexpidx));
    z_perline = accumarray(exp2lineidx(~badexpidx)',expweights(~badexpidx));
    muhat_perline = sum_perline ./ z_perline;
    muhat_perline = muhat_perline';
    
    % given the current values for muhat_perline, estimate sigmahat_set and
    % sigmahat_exp (keeping weights of sets constant)
    nsets = numel(setstats.metadata);
    sigmahat_set_totals = nan(1,nsets);
    sigmahat_exps = nan(1,nsets);
    for seti = 1:nsets,
      
      expidx_curr = find(exp2set == seti & ~badexpidx);
      linei = set2lineidx(seti);
      ncurr = nexpsperset(seti);
      if ncurr <= 2 || badsetidx(seti),
        continue;
      end
      mucurr = nanmean(x(expidx_curr));
      sigmahat_exps(seti) = sqrt(sum( (x(expidx_curr)-mucurr).^2 )/(numel(expidx_curr)-1));
      sigmahat_set_totals(seti) = (mucurr - muhat_perline(linei))^2;
      
    end
    idx = nexpsperset > 2;
    sigmahat_set_total = sqrt( sum(sigmahat_set_totals(idx)) / (nnz(idx)-1) );
    sigmahat_exp_new = sqrt(mean(sigmahat_exps(idx).^2));
    mean_nexpsperset = mean(nexpsperset(idx));
    if sigmahat_exp_new >= sigmahat_set_total,
      sigmahat_set_new = 0;
    else
      sigmahat_set_new = sqrt(sigmahat_set_total^2 - sigmahat_exp_new^2 / mean_nexpsperset);
    end
    
    maxdiff = max(abs(sigmahat_set_new-sigmahat_set),abs(sigmahat_exp_new-sigmahat_exp));
    
    sigmahat_set = sigmahat_set_new;
    sigmahat_exp = sigmahat_exp_new;
    
    fprintf('iter %d: sigmahat_set = %f, sigmahat_exp = %f\n',iter,sigmahat_set,sigmahat_exp);
    
    if maxdiff < .01,
      break;
    end
    
  end
  
  fprintf('%s; after %d iterations, chose sigmahat_set = %f and sigmahat_exp = %f\n',statfn,iter,sigmahat_set,sigmahat_exp);
  if sigmahat_exp / sigmahat_set > 100,
    sigmahat_set = sigmahat_exp / 100;
    fprintf('Ratio between sigmahat_exp and sigmahat_set too large, setting it to 100\n');
  end
  sigmahat_set_perstat(stati) = sigmahat_set;
  sigmahat_exp_perstat(stati) = sigmahat_exp;
  
  expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
  expweights = expweights / sum(expweights);
  
  sum_perline = accumarray(exp2lineidx(~badexpidx)',x(~badexpidx)'.*expweights(~badexpidx));
  sum2_perline = accumarray(exp2lineidx(~badexpidx)',x(~badexpidx).^2'.*expweights(~badexpidx));
  z_perline = accumarray(exp2lineidx(~badexpidx)',expweights(~badexpidx));
  muhat_perline = sum_perline ./ z_perline;
  sigmahat_perline = sqrt(sum2_perline ./ z_perline - muhat_perline.^2);
  
  muhat_weighted_perstat(:,stati) = muhat_perline;
  sigmahat_weighted_perstat(:,stati) = sigmahat_perline;
  
end

%% randomly construct similar-sized data sets

nexpsperset_control = setstats.nexps.(statfn)(setiscontrol);
set2exp_control = cell(1,numel(setidxcontrol));
for i = 1:numel(setidxcontrol),
  set2exp_control{i} = find(exp2set == setidxcontrol(i));
end

[~,order] = sort(cellfun(@numel,statfns(statis)));

nbiggersamples = nan(nlines,nstats);
nsmallersamples = nan(nlines,nstats);
nmoreextremesamples = nan(nlines,nstats);

for stati = statis(:)',
  
  statfn = statfns{stati};
  maxnexps = max(setstats.nexps.(statfn));
  
  x = allstats.(statfn);
  badexpidx = isnan(x);
  
  sigmahat_set = sigmahat_set_perstat(stati);
  sigmahat_exp = sigmahat_exp_perstat(stati);
  
  for linei = 1:nlines,
    
    nsets_curr = linestats.nsets.(statfn)(linei);
    setidx_curr = find(set2lineidx == linei);
    nexpsperset_curr = setstats.nexps.(statfn)(setidx_curr);
    absdiffline = abs(muhat_weighted_perstat(linei,stati) - muhat_weighted_perstat(linecontrolidx,stati));
    
    counts = hist(nexpsperset_curr,1:maxnexps);
    
    nbiggersamples(linei,stati) = 0;
    nsmallersamples(linei,stati) = 0;
    nmoreextremesamples(linei,stati) = 0;
    for samplei = 1:nsamples,
      
      setidxsample = [];
      expidxsample = [];
      isleft = true(1,numel(setidxcontrol));
      n = [];
      for minnexps = maxnexps:-1:1,
        
        % choose some sets
        isallowed = nexpsperset_control >= minnexps;
        idxsample = randsample(find(isleft&isallowed),counts(minnexps),false);
        isleft(idxsample) = false;
        
        setidxsample_curr = setidxcontrol(idxsample);
        
        % choose some experiments
        for seti = 1:counts(minnexps),
          expidxsample(end+1:end+minnexps) = randsample(set2exp_control{idxsample(seti)},minnexps);
        end
        
        n(end+1:end+counts(minnexps)*minnexps) = minnexps;
        
      end
      
      badidx = isnan(x(expidxsample));
      expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
      musample = sum(x(expidxsample(~badidx)).*expweights(~badidx)) / sum(expweights(~badidx));
      
      if musample < muhat_weighted_perstat(linei,stati),
        nsmallersamples(linei,stati) = nsmallersamples(linei,stati) + 1;
      elseif musample > muhat_weighted_perstat(linei,stati),
        nbiggersamples(linei,stati) = nbiggersamples(linei,stati) + 1;
      end
      
      absdiffsample = abs(musample - muhat_weighted_perstat(linecontrolidx,stati));
      if absdiffsample > absdiffline,
        nmoreextremesamples(linei,stati) = nmoreextremesamples(linei,stati) + 1;
      end
      
    end
    
    fprintf('%s, %s:\n',statfn,linestats.line_names{linei});
    fprintf('Percentile of %s mean among sample means: %.1f%%\n',...
      linestats.line_names{linei},nsmallersamples(linei,stati)/nsamples*100);
    fprintf('Percent of samples more different from control mean than %s: %.1f%%\n',...
      linestats.line_names{linei},nmoreextremesamples(linei,stati)/nsamples*100);
    
  end
  
  save pvaluescurrmb.mat nsmallersamples nbiggersamples nmoreextremesamples linei stati;
  
end

fracbiggersamples = nbiggersamples / nsamples;
fracsmallersamples = nsmallersamples / nsamples;
fracmoreextremesamples = nmoreextremesamples / nsamples;

%% get more precision if necessary

nsamples_dynamic = repmat(nsamples,[nlines,nstats]);
nsamples_curr = 10^4;

for stati = statis',
  
  statfn = statfns{stati};
  maxnexps = max(setstats.nexps.(statfn));
  
  x = allstats.(statfn);
  badexpidx = isnan(x);
  
  sigmahat_set = sigmahat_set_perstat(stati);
  sigmahat_exp = sigmahat_exp_perstat(stati);
    
  lineis = find(fracsmallersamples(:,stati)>=.999 | fracsmallersamples(:,stati)<=.001);
  if isempty(lineis),
    continue;
  end
  
  for linei = lineis(:)',
    
    % make sure it is possible to get more precision
    if fracsmallersamples(linei,stati) == 0,
      if min(x(expiscontrol)) > muhat_weighted_perstat(linei,stati),
        continue;
      end
    else
      if max(x(expiscontrol)) < muhat_weighted_perstat(linei,stati),
        continue;
      end
    end
    
    nsets_curr = linestats.nsets.(statfn)(linei);
    setidx_curr = find(set2lineidx == linei);
    nexpsperset_curr = setstats.nexps.(statfn)(setidx_curr);
    absdiffline = abs(muhat_weighted_perstat(linei,stati) - muhat_weighted_perstat(linecontrolidx,stati));
    
    counts = hist(nexpsperset_curr,1:maxnexps);
    
    nbiggersamples(linei,stati) = 0;
    nsmallersamples(linei,stati) = 0;
    nmoreextremesamples(linei,stati) = 0;
    nsamples_dynamic(linei,stati) = nsamples_curr;
    for samplei = 1:nsamples_curr,
      
      setidxsample = [];
      expidxsample = [];
      isleft = true(1,numel(setidxcontrol));
      n = [];
      for minnexps = maxnexps:-1:1,
        
        % choose some sets
        isallowed = nexpsperset_control >= minnexps;
        idxsample = randsample(find(isleft&isallowed),counts(minnexps),false);
        isleft(idxsample) = false;
        
        setidxsample_curr = setidxcontrol(idxsample);
        
        % choose some experiments
        for seti = 1:counts(minnexps),
          expidxsample(end+1:end+minnexps) = randsample(set2exp_control{idxsample(seti)},minnexps);
        end
        
        n(end+1:end+counts(minnexps)*minnexps) = minnexps;
        
      end
      
      badidx = isnan(x(expidxsample));
      expweights = 1./(n.*sigmahat_set^2 + sigmahat_exp^2);
      musample = sum(x(expidxsample(~badidx)).*expweights(~badidx)) / sum(expweights(~badidx));
      
      if musample < muhat_weighted_perstat(linei,stati),
        nsmallersamples(linei,stati) = nsmallersamples(linei,stati) + 1;
      elseif musample > muhat_weighted_perstat(linei,stati),
        nbiggersamples(linei,stati) = nbiggersamples(linei,stati) + 1;
      end
      
      absdiffsample = abs(musample - muhat_weighted_perstat(linecontrolidx,stati));
      if absdiffsample > absdiffline,
        nmoreextremesamples(linei,stati) = nmoreextremesamples(linei,stati) + 1;
      end
      
    end
    
    fprintf('%s, %s:\n',statfn,linestats.line_names{linei});
    fprintf('Percentile of %s mean among sample means: %.1f%%\n',...
      linestats.line_names{linei},nsmallersamples(linei,stati)/nsamples_curr*100);
    fprintf('Percent of samples more different from control mean than %s: %.1f%%\n',...
      linestats.line_names{linei},nmoreextremesamples(linei,stati)/nsamples_curr*100);
    
  end
  
  save pvaluescurrmb.mat nsmallersamples nbiggersamples nmoreextremesamples nsamples_dynamic linei stati;
  
end

fracsmallersamples = nsmallersamples ./ nsamples_dynamic;
fracbiggersamples = nbiggersamples ./ nsamples_dynamic;

%% choose the most significant test, and correct for this

pvalue = min(1,min(1-fracsmallersamples,1-fracbiggersamples)*2);
min_pvalue_observable = 2 / max(nsamples_dynamic(:));

%% make a table

statfns_plot = {
  'fractime_flyany_framestop'
  'fractime_flyany_framewalk'
  'fractime_flyany_framejump'
  'fractime_flyany_framerighting'
  'fractime_flyany_framepivottail'
  'fractime_flyany_framepivotcenter'
  'fractime_flyany_framebodyturn'
  'fractime_flyany_framebackup'
  'fractime_flyany_framecrabwalkall'
  'fractime_flyany_framecrabwalkextreme'
  'fractime_flyany_framechase'
  'fractime_flyany_frameattemptedcopulation'};

% sort by mean for statfns_plot{1};
stati = find(strcmp(statfns,statfns_plot{1}));
[~,lineorder] = sort(muhat_weighted_perstat(:,stati),'descend');

% normalize strength
ismaybesig = pvalue < .05;
dcontrol = bsxfun(@minus,muhat_weighted_perstat,muhat_weighted_perstat(linecontrolidx,:));
tmp = dcontrol;
tmp(~ismaybesig) = nan;
minv = min(tmp,[],1);
maxv = max(tmp,[],1);

log10pvalue = log10(max(min_pvalue_observable,pvalue));

ncolors = 256;
cm = redblue(2*ncolors);
% cm = zeros(2*ncolors,3);
% cm(:,3) = linspace(1,0,2*ncolors);
% cm(:,1) = linspace(0,1,2*ncolors);
% cm(1:ncolors,2) = linspace(0,1,ncolors);
% cm(ncolors+1:end,2) = linspace(1,0,ncolors);
cm = bsxfun(@rdivide,cm,sqrt(sum(cm.^2,2)));

isneg = dcontrol < 0;
coloridx = nan(size(dcontrol));
z = -minv;
z(minv==0) = 1;
tmp = max(1,floor(bsxfun(@rdivide,bsxfun(@minus,dcontrol,minv),z)*ncolors)+1);
coloridx(isneg) = tmp(isneg);
z = maxv;
z(maxv==0) = 1;
tmp = min(ncolors,floor(bsxfun(@rdivide,dcontrol,z)*ncolors)+1)+ncolors;
coloridx(~isneg) = tmp(~isneg);
coloridx = min(2*ncolors,max(1,coloridx));

hue = reshape(cm(coloridx,:),[nlines,nstats,3]);

maxp = max(log10pvalue(:));
minp = min(log10pvalue(:));
intensity = 1-(log10pvalue-minp)/(maxp-minp);
color = 1-bsxfun(@times,1-hue,intensity);

color_ordered = color(lineorder,:,:);

[ism,idx] = ismember(statfns,statfns_plot);
statis = find(ism);
[~,order] = sort(idx(ism));
statis = statis(order);

color_ordered = color_ordered(:,statis,:);

hfig = 1;
figure(hfig);
clf;
image(permute(color_ordered,[2,1,3]));

behaviornames = regexprep(statfns_plot,'^fractime_flyany_frame','');
set(gca,'YTick',1:numel(statis),'YTickLabels',behaviornames);

short_line_names = regexprep(linestats.line_names,'^GMR_','');
short_line_names = regexprep(short_line_names,'_AE_01$','');
short_line_names = regexprep(short_line_names,'_01$','');
short_line_names = regexprep(short_line_names,'_BX3G_1500154$','');

set(gca,'XTick',1:nlines,'XTickLabel',short_line_names(lineorder),'box','off');