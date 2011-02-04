% choose the area thresholds for male vs female

%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;

%% data locations

protocol = '20110111';
settingsdir = 'settings';
dataloc_paramsfilestr = 'dataloc_params.txt';

[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol',['CtraxTest',protocol]);

dataloc_params = ReadParams(fullfile(settingsdir,protocol,dataloc_paramsfilestr));

sexclassifiermatfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/sexclassifier.mat';

%% load the data

trx = Trx('rootreaddir',dataloc_params.rootreaddir,...
  'rootwritedir',dataloc_params.rootwritedir,...
  'perframedir',dataloc_params.perframedir,...
  'trxfilestr',dataloc_params.trxfilestr,...
  'sexclassifiermatfile',sexclassifiermatfile);

for i = 1:numel(expdirs),
  moviefile = fullfile(rootreaddir,expdirs{i},dataloc_params.moviefilestr);
  [~,~,fid,vidinfo] = get_readframe_fcn(moviefile);
  fclose(fid);
  trx.AddExpDir(expdirs{i},vidinfo);
end

%% filter out areas that are maybe caused by tracking noise

maxfreq = .005;
filterorder = 1;
outlierprctile = .1;
mindur = 20;

f = fdesign.lowpass('N,F3db',filterorder,maxfreq);
h = design(f,'butter');
h.PersistentMemory = true;

errx = [];
areasmooth = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  h.filter(fliplr(trx(fly).area));
  areasmooth{fly} = h.filter(trx(fly).area);
  errx = [errx,abs(areasmooth{fly}-trx(fly).area)]; %#ok<AGROW>
end

maxerrx = prctile(errx,100-outlierprctile);

isoutlierarea = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  isoutlierarea{fly} = abs(areasmooth{fly}-trx(fly).area) > maxerrx;
end


%% learn a 2-state HMM for area in an unsupervised manner

nstates = 2;
state2sex = cell(1,nstates);
if true,
  % break into smaller sequences where we're sure of area estimate
  X = cell(1,trx.nflies);
  for fly = 1:trx.nflies,
    areacurr = trx(fly).area;
    % interpolate across bad data
    [starts,ends] = get_interval_ends(isoutlierarea{fly});
    ends = ends - 1;
    for i = 1:numel(starts),
      if starts(i) == 1,
        areacurr(starts(i):ends(i)) = areacurr(ends(i)+1);
      elseif ends(i) == trx(fly).nframes,
        areacurr(starts(i):ends(i)) = areacurr(starts(i)-1);
      else
        areacurr(starts(i):ends(i)) = (areacurr(starts(i)-1)+areacurr(ends(i)+1))/2;
      end
    end
    X{fly} = areacurr(:);
  end
%   X = trx.area;
%   for i = 1:numel(X),
%     X{i} = X{i}';
%   end
  [mu_area,var_area,ptrans,prior,ll]=hmm_multiseq(X,nstates);
  state2sex{argmax(mu_area)} = 'F';
  state2sex{argmin(mu_area)} = 'M';
  save sexclassifier.mat mu_area var_area ptrans prior ll nstates state2sex maxerrx;
else
  load sexclassifier.mat;
end

%% apply the hmm

x_trans_fun = @(areacurr,areaprev,scurr) normpdf(areacurr,mu_area(scurr),var_area);
log_ptrans = log(ptrans);
log_s_trans_fun = @(scurr,sprev) log_ptrans(sprev,scurr);
for fly = 1:trx.nflies,
  sexstate = first_order_viterbi(trx(fly).areasmooth',nstates,x_trans_fun,log_s_trans_fun);
  trx.SetSex(fly,state2sex(sexstate));
end

%%