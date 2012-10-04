% ScriptControlNormalizationTest20120727

load AllExperimentsFracTimePerSet20120601T141100.mat;
nbehaviors = size(fractime,2);

%%

[sets,~,setidx] = unique({metadata.set});
muset = nan(numel(sets),nbehaviors);
sigmaset = nan(numel(sets),nbehaviors);
nset = nan(1,numel(sets));
for i = 1:numel(sets),
  muset(i,:) = mean(fractime(setidx==i,:),1);
  sigmaset(i,:) = std(fractime(setidx==i,:),1,1);
  nset(i) = nnz(setidx==i);
end

idxcontrol = ~cellfun(@isempty,regexp(sets,'^pBD','once'));

%%

set_datetimes = regexp(sets,'\d{8}T\d{6}$','once','match');
set_datenums = datenum(set_datetimes,'yyyymmddTHHMMSS')';
set_day = floor(set_datenums);
days = unique(set_day);

muday = nan(numel(days),nbehaviors);
sigmaday = nan(numel(days),nbehaviors);
for i = 1:numel(days),
  idxcurr = days(i) == set_day & idxcontrol;
  if ~any(idxcurr), continue; end
  muday(i,:) = mean(muset(idxcurr,:),1);
  sigmaday(i,:) = std(muset(idxcurr,:),1,1);
end

gmrmuday = nan(numel(days),nbehaviors);
gmrsigmaday = nan(numel(days),nbehaviors);
for i = 1:numel(days),
  idxcurr = days(i) == set_day & ~idxcontrol;
  if ~any(idxcurr), continue; end
  gmrmuday(i,:) = mean(muset(idxcurr,:),1);
  gmrsigmaday(i,:) = std(muset(idxcurr,:),1,1);
end


%% plot per-day means for gmrs and controls

hfig = 1;
figure(hfig);
clf;
hax = createsubplots(1,nbehaviors,.05);

for behi = 1:nbehaviors,
  axes(hax(behi));
  plot(repmat(days+.15,[2,1]),bsxfun(@plus,gmrmuday(:,behi),bsxfun(@times,[-1,1],gmrsigmaday(:,behi)))','-','color',[.5,.5,.75]);  
  hold on;
  plot(repmat(days-.15,[2,1]),bsxfun(@plus,muday(:,behi),bsxfun(@times,[-1,1],sigmaday(:,behi)))','-','color',[.75,.5,.5]);
  hgmr = plot(set_day(~idxcontrol)+.15,muset(~idxcontrol,behi),'.','color',[0,0,.7]);
  hcontrol = plot(set_day(idxcontrol)-.15,muset(idxcontrol,behi),'.','color',[.7,0,0]);
  set(hax(behi),'XLim',[days(1)-1,days(end)+1]);
  ylim = get(hax(behi),'YLim');
  ylim(1) = -diff(ylim)/40;
  set(hax(behi),'YLim',ylim);
  datetick('x','yymmmdd','keeplimits','keepticks');
  title(scoresfilestr{behi,1});
end

legend([hgmr,hcontrol],{'GMR','Control'});

%% per-day means for controls vs gmrs

hfig = 2;
figure(hfig);
clf;
hax = createsubplots(1,nbehaviors,.05);

for behi = 1:nbehaviors,
  axes(hax(behi));
  plot(muday(:,behi),gmrmuday(:,behi),'k.');
  title(scoresfilestr{behi,1});
end

%% 

set_line_names = regexp(sets,'^(.*)__Rig','tokens','once');
set_line_names = cellfun(@(x) x{1},set_line_names,'UniformOutput',false);

[unique_line_names,~,lineidx] = unique(set_line_names);
idxignore = ismember(unique_line_names,{'pBDPGAL4U','FCF_pBDPGAL4U_1500437'});

nsetsperline = hist(lineidx,1:numel(unique_line_names));
nsetsperline(idxignore) = 0;

%% 

radpools = 0:30;

for pooli = 1:numel(radpools),
  radpool = radpools(pooli);

  fil = ones(2*pooli+1,1)/(2*pooli+1);
  % TODO: if we do this for real, we want to be more careful at the boundaries
  controlnorm = imfilter(muday,'replicate');
  

muline = nan(numel(unique_line_names),nbehaviors);
sigmaline = nan(numel(unique_line_names),nbehaviors);

for i = 1:numel(unique_line_names),
  if idxignore(i),
    continue;
  end
  muline(i,:) = mean(muset(lineidx==i,:),1);
  sigmaline(i,:) = std(muset(lineidx==i,:),1,1);
  nsetsperline(i) = nnz(lineidx==i);
end

