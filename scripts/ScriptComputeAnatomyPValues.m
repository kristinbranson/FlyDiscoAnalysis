%% set up paths

anatomydatadir = '/nobackup/branson/AverageAnatomyData20130618';
metadatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryMetadata20130912.mat';
maskfile = '/groups/branson/bransonlab/projects/olympiad/cross_assay/trunk/matlab/kristin/FullBrainMask.mat';

%% parameters

load(maskfile);
mask0 = permute(mask,[2,1,3]);
mask = mask0 > 0;

lowthresh = .4;
highthresh = .6;
minnlineswithdata = 500;
alphabinofit = .001;

%% load metadata

load(metadatafile);
line_names_behavior = unique({metadata.line_name});

tmp = dir(fullfile(anatomydatadir,'meanim*.mat'));
line_names_anatomy = cellfun(@(x) x{1},regexp({tmp.name},'^meanim_(.*)\.mat$','tokens','once'),'UniformOutput',false);

line_names = intersect(line_names_behavior,line_names_anatomy);
nlines = numel(line_names);

%% 

hfig = 1;
figure(hfig);
clf;
hax = gca;

i = 1;
filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
datacurr = matfile(filename);
[nr,nc,nz] = size(datacurr.meanim);
npx = nr*nc*nz;
countspos = zeros([1,npx]);
countsneg = zeros([1,npx]);

lineorder = randperm(nlines);
tmpedges = linspace(0,1,101);
tmpcenters = (tmpedges(1:end-1)+tmpedges(2:end))/2;

for ii = 1:nlines,
  fprintf('Line %d / %d\n',ii,nlines);
  i = lineorder(ii);
  filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
  datacurr = matfile(filename);
  idx = (datacurr.meanim>highthresh)&mask;
  countspos(idx) = countspos(idx) + 1;
  idx = (datacurr.meanim<lowthresh)&mask;
  countsneg(idx) = countsneg(idx) + 1;
  
  if mod(ii,100) == 0,
    
    if ~ishandle(hax),
      figure(hfig);
      clf;
      hax = gca;
    end
    hold(hax,'off');
    tmppos = hist(countspos/max(countspos(:)),tmpcenters);
    tmpneg = hist(countsneg/max(countsneg(:)),tmpcenters);
    
    plot(hax,tmpcenters,tmpneg,'k.-');
    hold(hax,'on');
    plot(hax,tmpcenters,tmppos,'r.-');
    set(hax,'YScale','log');
    title(num2str(ii));
    drawnow;
    
    save tmp.mat countsneg countspos highthresh lowthresh line_names lineorder;
  end
end

nlinespervoxel = reshape(countspos+countsneg,[nr,nc,nz]);
fraclineswithexpression = reshape(countspos+1,[nr,nc,nz])./(nlinespervoxel+1);
isdata = nlinespervoxel >= minnlineswithdata & mask;

% get upper and lower bounds on the bias
pmean = zeros([nr,nc,nz]);
plowerbounds = zeros([nr,nc,nz]);
pupperbounds = zeros([nr,nc,nz]);
for z = 1:nz,
  fprintf('%d / %d\n',z,nz);
  [tmp1,tmp2] = binofit(reshape(countspos(:,:,z),[nr*nc,1]),reshape(countspos(:,:,z)+countsneg(:,:,z),[nr*nc,1]),alphabinofit);
  pmean(:,:,z) = reshape(tmp1,[nr,nc]);
  plowerbounds(:,:,z) = reshape(tmp2(:,1),[nr,nc]);
  pupperbounds(:,:,z) = reshape(tmp2(:,2),[nr,nc]);
end


save -v7.3 AllLinesAnatomyPValueData20130917.mat fraclineswithexpression pmean plowerbounds pupperbounds nlinespervoxel isdata countspos countsneg lowthresh highthresh line_names;

%% compute a table of pvalues because computing p-values takes a long time!

ps = linspace(min(plowerbounds(:)),max(pupperbounds(:)),100);
ns = 1:nlines;
xs = 0:nlines;
minp = ps(1);
maxp = ps(end);
dp = mean(diff(ps));

pvaluetable = zeros([numel(xs),numel(ns),numel(ps)]);
for pi = 1:numel(ps),
  fprintf('%d/%d...\n',pi,numel(ps));
  for ni = 1:numel(ns),
    idx = xs <= ns(ni) & xs > 0;
    pvaluetable(idx,ni,pi) = binocdf(xs(idx)-1,ns(ni),ps(pi));
  end
end
pvaluetable = 1-pvaluetable;

nps = numel(ps);
nns = numel(ns);
nxs = numel(xs);
a = (nps-1)/(maxp-minp);
pvaluetablefun = @(x,n,p) pvaluetable(sub2ind([nxs,nns,nps],x+1,n,max(1,ceil( 1+(p-minp)*a ))));

save -append AllLinesAnatomyPValueData20130917.mat ps ns xs minp maxp dp pvaluetable nps nns nxs a pvaluetablefun;

%% compute p-values for random sets of lines

i = 1;
filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
datacurr = matfile(filename);
[nr,nc,nz] = size(datacurr.meanim);

nsample = 10;
linessample = randsample(nlines,nsample)';
countsposcurr = zeros([nr,nc,nz]);
countsnegcurr = zeros([nr,nc,nz]);
for ii = 1:nsample,
  fprintf('Sample %d / %d\n',ii,nsample);
  i = linessample(ii);
  filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
  datacurr = matfile(filename);
  countsposcurr(datacurr.meanim>highthresh) = countsposcurr(datacurr.meanim>highthresh) + 1;
  countsnegcurr(datacurr.meanim<lowthresh) = countsnegcurr(datacurr.meanim<lowthresh) + 1;
end

pvaluecurr = pvaluetablefun(countsposcurr,countsposcurr+countsnegcurr,pupperbounds);
pvaluecurr(~mask) = 1;

%% estimate false positive rates for various p-values, nsamples

[nr,nc,nz] = size(pupperbounds);

edges = [0,10.^(-16:-2),.02,.03,.04,.05,1];
centers = (edges(2:end)+edges(1:end-1))/2;

maxn = 100;
nsamples = 10;

countsposcurr = zeros([nr,nc,nz]);
countsnegcurr = zeros([nr,nc,nz]);
pvaluecurr = ones([nr,nc,nz]);

frac = zeros(numel(edges),maxn,nsamples);

colors = (jet(maxn)+.25)/1.25;
hfig = 1;
figure(hfig);
clf;
hax = gca;
set(hax,'Color','k','XScale','log');
hold on;

for samplei = 1:nsamples,

  countsposcurr(:) = 0;
  countsnegcurr(:) = 0;
  linessample = randsample(nlines,maxn)';

  for n = 1:maxn

    fprintf('Reading in image %d / %d for sample %d / %d\n',n,maxn,samplei,nsamples);
    i = linessample(n);
    filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_names{i}));
    datacurr = load(filename);
    countsposcurr(datacurr.meanim>highthresh) = countsposcurr(datacurr.meanim>highthresh) + 1;
    countsnegcurr(datacurr.meanim<lowthresh) = countsnegcurr(datacurr.meanim<lowthresh) + 1;
    
    goodidx = (countsposcurr+countsnegcurr)>0;

    pvaluecurr(:) = 1;
    maskcurr = goodidx&mask;
    pvaluecurr(maskcurr) = pvaluetablefun(countsposcurr(maskcurr),countsposcurr(maskcurr)+countsnegcurr(maskcurr),pupperbounds(maskcurr));
    frac(:,n,samplei) = histc(pvaluecurr(maskcurr),edges) / nnz(maskcurr);
    
    plot(hax,[centers,1],frac(:,n,samplei),'.-','Color',colors(n,:));
    set(hax,'YLim',[-.01,max(frac(:))*1.01]);
    drawnow;

  end
end