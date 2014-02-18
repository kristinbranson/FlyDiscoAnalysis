
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStats20131127.mat';
outmatfile = 'ComputePValueBySamplingData20131127.mat';
resultsdir = 'ComputePValueBySamplingResults20131127';
  
%%  load in data

load(datafile);

%% compare to randomly sampled control data

% which sets are control
setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
setidxcontrol = find(setidxcontrol);

% experiment to set index
[~,exp2setidx] = ismember({metadata.set},{setstats.metadata.set});

% experiments that are control
expiscontrol = find(ismember(exp2setidx,setiscontrol));

% number of control sets
ncontrolsets = nnz(setiscontrol);

% number of experiments in each control set
statfn = 'velmag_flyany_frameany';
setnexps_control = setstats.nexps.(statfn)(setiscontrol)';

controldatenum = floor(datenum({setstats.metadata(setiscontrol).exp_datetime},'yyyymmddTHHMMSS'));

nsamples = 100000;

maxnexps = max(setstats.nexps.(statfn));

minnframes = repmat(200,[1,nstats]);
minnexps = 2;

save(outmatfile,...
  'statfns','nstats','setstats','allstats','allnframestotal','minnframes','minnexps',...
  'setiscontrol','nsets','maxnexps','exp2setidx','nlines','linestats','set2lineidx',...
  'nsamples');


%% run jobs on the cluster
%

fracsmaller = nan(nlines,nstats);
fracbigger = nan(nlines,nstats);
for stati = 1:nstats,
  [fracsmaller(1:end-1,stati),fracbigger(1:end-1,stati)] = ComputePValueBySampling(stati,'ComputePValueBySamplingData20130925.mat','ComputePValueBySamplingResults20130925');
end
% 
% 
% for stati = stati:nstats,
%   
%   statfn = statfns{stati};
%   
%   fprintf('Statistic %s: %d / %d\n',statfn,stati,nstats);
%   
%   isgoodset = ~isinf(setstats.normmeans.(statfn)) & ...
%     ~isnan(setstats.normmeans.(statfn)) & ...
%     setstats.nexps.(statfn) >= minnexps;
%   isgoodexp = ~isinf(allstats.(statfn)) & ...
%     ~isnan(allstats.(statfn)) & ...
%     allnframestotal.(statfn) >= minnframes(stati);
%   
%   isgoodcontrolset = isgoodset & setiscontrol;
%   if nnz(isgoodcontrolset <= 1),
%     fprintf('One or less good control sets for statistic %s, not computing p-values\n',statfn);
%     continue;
%   end
%   
%   
%   idxgoodcontrolset = find(isgoodcontrolset);
%   
%   set2expcurr = nan(nsets,maxnexps);
%   for seti = idxgoodcontrolset,
%     tmp = find(exp2setidx==seti & isgoodexp);
%     set2expcurr(seti,1:numel(tmp)) = tmp;
%   end
% 
%   parfor linei = 1:nlines,
% 
%     if mod(linei,100) == 0,
%       fprintf('Stat %s (%d / %d), line %s (%d / %d)\n',statfn,stati,nstats,...
%         linestats.line_names{linei},linei,nlines);
%     end
%     
%     %fprintf('Line %s: %d / %d\n',linestats.line_names{linei},linei,nlines);
%   
%     % number of experiments in each set
%     setidxcurr = find(set2lineidx==linei & isgoodset);
%     if isempty(setidxcurr),
%       %fprintf('No good sets found for line %s, skipping.\n',linestats.line_names{linei});
%       continue;
%     end
%     
%     nexpscurr = setstats.nexps.(statfn)(setidxcurr);
%     datenumscurr = setdatenum(set2lineidx==linei & isgoodset);
%     nexpscurr = sort(nexpscurr,'descend');
%     nsetscurr = numel(nexpscurr);
%     
%     % sample
%     setnormmu = nan(nsamples,nsetscurr);
%     for setii = 1:nsetscurr,
%       
%       setidxallowed = find(isgoodcontrolset & ...
%         setstats.nexps.(statfn) >= nexpscurr(setii));
%       for nsub = 1:nexpscurr(setii)-1,
%         if numel(setidxallowed) > 1,
%           break;
%         end
%         setidxallowed = find(isgoodcontrolset & ...
%           setstats.nexps.(statfn) >= nexpscurr(setii)-nsub);
%       end
%       nsub = nsub-1;
%       if nsub > 0,
%         fprintf('For stat %s, line %s, set %d, needed to consider sets %d exps smaller than than this set\n',statfn,linestats.line_names{linei},setii,nsub);
%         nexpscurr(setii) = nexpscurr(setii)-nsub;
%       end
%       setsampleis = randsample(setidxallowed,nsamples,true);
%             
%       % choose the experiments per set
%       controlnexpscurr = setstats.nexps.(statfn)(setsampleis);
%       
%       % by default, use the first experiments
%       expsampleis = set2expcurr(setsampleis,1:nexpscurr(setii));
%       
%       % for sets with more experiments, sample without replacement
%       tmpidx = find(controlnexpscurr > nexpscurr(setii));
%       for tmpi = tmpidx,
%         expsampleis(tmpi,:) = set2expcurr(setsampleis(tmpi),randsample(controlnexpscurr(tmpi),nexpscurr(setii)));
%       end
%       
%       tmp = allstats.(statfn)(expsampleis);
%       tmp = reshape(tmp,[nsamples,nexpscurr(setii)]);
%       mu = mean(tmp,2);
%       setnormmu(:,setii) = mu - setstats.controlmeans.(statfn)(setsampleis)';
% 
%     end
% 
%     mu = mean(setnormmu,2);
%     fracsmaller(linei,stati) = nnz(mu<linestats.normmeans.(statfn)(linei))/nsamples;
%     fracbigger(linei,stati) = nnz(mu>linestats.normmeans.(statfn)(linei))/nsamples;
%   
%   end
%   
% end

%% load results run on cluster

fracsmaller = nan(nlines,nstats);
fracbigger = nan(nlines,nstats);
for stati = 1:nstats,
  
  filename = fullfile(resultsdir,sprintf('PvaluesForStat%03d.mat',stati));
  if ~exist(filename,'file'),
    fprintf('Stat %d failed\n',stati);
    continue;
  end
  
  tmp = load(filename);
  fracsmaller(1:end-1,stati) = tmp.fracsmaller_stat;
  fracbigger(1:end-1,stati) = tmp.fracbigger_stat;
  
end

% remove p-values for missing data
ismissing = false(nlines,nstats);
for stati = 1:nstats,
  ismissing(:,stati) = isnan(linestats.normmeans.(statfns{stati})) | ...
    isinf(linestats.normmeans.(statfns{stati}));
end

fracsmaller(ismissing) = nan;
fracbigger(ismissing) = nan;


save -append CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat fracsmaller fracbigger

%% compute standard deviation of control line for each statistic
% controllinestd = nan(1,nstats);
% controllinemean = nan(1,nstats);
% for stati = 1:nstats,
%   ComputeControlLineStdBySampling;
%   controllinemean(stati) = mucurr;
%   controllinestd(stati) = sigcurr;
% end

%% compute upper and lower bounds on these p-values at 99.9% confidence

alpha = .0001;
[tmp1,tmp2] = binofit((1-fracbigger(:))*nsamples,nsamples,alpha);
tmp1 = reshape(tmp1,[nlines,nstats]);
fracsmaller_upperbound = reshape(tmp2(:,2),[nlines,nstats]);
pvalue_smaller = min(1,fracsmaller_upperbound+alpha);
pvalue_smaller(ismissing) = nan;

[tmp1,tmp2] = binofit((1-fracsmaller(:))*nsamples,nsamples,alpha);
tmp1 = reshape(tmp1,[nlines,nstats]);
fracbigger_upperbound = reshape(tmp2(:,2),[nlines,nstats]);
pvalue_bigger = min(1,fracbigger_upperbound+alpha);
pvalue_bigger(ismissing) = nan;

save -append CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat pvalue_smaller pvalue_bigger fracsmaller_upperbound fracbigger_upperbound;
