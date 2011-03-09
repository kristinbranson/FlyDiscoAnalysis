% choose the area thresholds for male vs female

%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath E:\Code\hmm
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
end

%% data locations

protocol = '20110202';
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
dataloc_paramsfilestr = 'dataloc_params.txt';

[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol',['CtraxTest',protocol],...
  'subwritefiles',{'registered_trx.mat'});

dataloc_params = ReadParams(fullfile(settingsdir,protocol,dataloc_paramsfilestr));

sexclassifierparamsfile = fullfile(settingsdir,protocol,dataloc_params.sexclassifierparamsfilestr);
sexclassifiermatfile = fullfile(settingsdir,protocol,dataloc_params.sexclassifiermatfilestr);
sexclassifiertxtfile = fullfile(settingsdir,protocol,dataloc_params.sexclassifiertxtfilestr);

sexclassifier_params = ReadParams(sexclassifierparamsfile);

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

maxfreq = sexclassifier_params.areasmooth_maxfreq;
filterorder = sexclassifier_params.areasmooth_filterorder;
outlierprctile = sexclassifier_params.areasmooth_outlierprctile;
mindur = sexclassifier_params.mindur;

f = fdesign.lowpass('N,F3db',filterorder,maxfreq);
h = design(f,'butter');
h.PersistentMemory = true;

errx = [];
areasmooth = cell(1,trx.nflies);
for fly = 1:trx.nflies,
  areasmooth{fly} = LowPassFilterArea(trx(fly).area,filterorder,maxfreq);
  errx = [errx,abs(areasmooth{fly}-trx(fly).area)]; %#ok<AGROW>
end

maxerrx = prctile(errx,100-outlierprctile);

%% learn a 2-state HMM for area in an unsupervised manner

nstates = 2;
state2sex = cell(1,nstates);
if true,
  % break into smaller sequences where we're sure of area estimate
  X = cell(1,trx.nflies);
  for fly = 1:trx.nflies,
    areacurr = SmoothAreaOutliers(trx(fly).area,filterorder,maxfreq,maxerrx);
    X{fly} = areacurr(:);
  end
  [mu_area,var_area,ptrans,prior,ll]=hmm_multiseq(X,nstates);
  state2sex{argmax(mu_area)} = 'F';
  state2sex{argmin(mu_area)} = 'M';
  save(sexclassifiermatfile,'mu_area','var_area','ptrans','prior','ll','nstates','state2sex','maxerrx',...
    'maxfreq','filterorder');
  
  fid = fopen(sexclassifiertxtfile,'w');
  ifemale = find(strcmp(state2sex,'F'));
  imale = find(strcmp(state2sex,'M'));
  fprintf(fid,'mu_area_female,%f\n',mu_area(ifemale));
  fprintf(fid,'mu_area_male,%f\n',mu_area(imale));
  fprintf(fid,'var_area,%f\n',var_area);
  fprintf(fid,'ptrans_female_given_female,%f\n',ptrans(ifemale,ifemale));
  fprintf(fid,'ptrans_male_given_female,%f\n',ptrans(ifemale,imale));
  fprintf(fid,'ptrans_female_given_male,%f\n',ptrans(imale,ifemale));
  fprintf(fid,'ptrans_male_given_male,%f\n',ptrans(imale,imale));
  fprintf(fid,'prior_female,%f\n',prior(ifemale));
  fprintf(fid,'prior_male,%f\n',prior(imale));
  fprintf(fid,'areasmooth_maxerrx,%f\n',maxerrx);
  fprintf(fid,'areasmooth_maxfreq,%f\n',maxfreq);
  fprintf(fid,'areasmooth_filterorder,%f\n',filterorder);
  fclose(fid);
  
else
  load(sexclassifiermatfile);
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