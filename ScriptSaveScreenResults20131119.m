%% parameters

load CollectedPrimaryPerFrameStats20131122.mat linestats setstats expstats;

main_control_line_name = 'pBDPGAL4U';

%%

fns = fieldnames(linestats.means);
savestuff = struct;
savestuff.expstats = expstats;
savestuff.setstats = rmfield(setstats,{'controlmeans','nsetscontrolnorm'});
savestuff.linestats = linestats;
%savestuff.linestats = rmfield(linestats,{'int_auto','dist_auto','int_manual','dist_manual'});

%% z-score

savestuff.setstats.zscores = struct;
savestuff.linestats.zscores = struct;
savestuff.expstats.zscores = struct;

for i = 1:numel(fns),
  fn = fns{i};
  controlsetdata = setstats.normmeans.(fn)(setiscontrol);
  mucontrol = nanmean(controlsetdata,2);
  sigcontrol = nanstd(controlsetdata,1,2);
  if sigcontrol <= 0,
    sigcontrol = 1;
  end
  savestuff.setstats.zscores.(fn) = (savestuff.setstats.normmeans.(fn) - mucontrol)/sigcontrol;
  savestuff.linestats.zscores.(fn) = (savestuff.linestats.normmeans.(fn) - mucontrol)/sigcontrol;
  savestuff.expstats.zscores.(fn) = (savestuff.expstats.normmeans.(fn) - mucontrol)/sigcontrol;

end

%% pvalues

savestuff.linestats.pvalue_bigger = struct;
savestuff.linestats.pvalue_smaller = struct;
for i = 1:numel(statfns),
  
  fn = statfns{i};
  savestuff.linestats.pvalue_bigger.(fn) = pvalue_bigger(:,i)';
  savestuff.linestats.pvalue_smaller.(fn) = pvalue_smaller(:,i)';
  
end

%% compute the standard deviations over experiments for the line
[~,exp2lineidx] = ismember({savestuff.expstats.metadata.line_name},{savestuff.linestats.metadata.line_name});

for i = 1:numel(statfns),
  
  fn = statfns{i};
  x = expstats.normmeans.(fn);
  
  idxcurr = ~isnan(x) & ~isinf(x) & exp2lineidx > 0;
  
  
  savestuff.linestats.normstds_over_exps.(fn) = accumarray(exp2lineidx(idxcurr)',x(idxcurr)',[nlines,1],@(y) std(y,1))';

  savestuff.linestats.pvalue_smaller.(fn) = pvalue_smaller(:,i)';
  
end


%%



timestamp = datestr(now,'yyyymmdd');

save(fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlResults',sprintf('FlyBowlScreenResults%s.mat',timestamp)),'-struct','savestuff');