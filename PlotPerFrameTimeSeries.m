function [hfig] = PlotPerFrameTimeSeries(trx,ind,field,varargin)

defaultcolors = struct('stimbkgd',[0.927 0.8156 0.8368],...
  'stim',[0.6350 0.0780 0.1840],...
  'off',[0 0 0],...
  'stimstd',[0.8175 0.539 0.592],...
  'offstd',[.7,.7,.7]...
);

% parse plotting parameters TODO add colors 
[hfig,visible,position,axposition,fontsize,bin,percbin,colors] = ...
  myparse(varargin,'hfig',[],...
  'visible','on',...
  'position',[1 1 800 400],...
  'axposition',[.075,.25,.9,.7],...
  'fontsize',8,...
  'bin',10,...
  'percbin',50,...
  'colors',defaultcolors);

fns = setdiff(fieldnames(defaultcolors),fieldnames(colors));
for i = 1:numel(fns),
  colors.(fns{i}) = defaultcolors.(fns{i});
end

% if percbah, load score data
ispercbeh = false;
if contains(field,'scores')
    ispercbeh = true;
    [~,data] = LoadScoresFromFile(trx,field,1);
end

% compute time series average and std for feature
%function [mean,std] = CopmuteTimeSeriesData(data)
nflies = trx.nflies;
% better way to get data size? 
nframes = max(trx.endframes);
timeseriesdata = nan(nflies,nframes);
%loop over fly data to add to matrix
dataflycurr =[];
for fly = 1:nflies
    if ispercbeh
        dataflycurr = data{fly};
    else
        dataflycurr = trx.GetPerFrameData(field,fly);
    end
    t0 = trx.firstframes(fly);
    t1 = t0 + numel(dataflycurr)-1;
    timeseriesdata(fly,t0:t1) = dataflycurr;
end


stimeseriesdata = nan(size(timeseriesdata));
if ispercbeh
    % compute percentage of flies behaving
    nbeh= sum(timeseriesdata,1,'omitnan');
    ntot = sum(~isnan(timeseriesdata),1);
    datamean = (nbeh./ntot).*100;
    % take moving average
    datamean= smooth(datamean',percbin,'moving')';
    datastd = zeros(1,numel(datamean));
else
    for fly = 1:nflies
        stimeseriesdata(fly,:) = smooth(timeseriesdata(fly,:)',bin,'moving')';
    end
    datamean = mean(stimeseriesdata,1,"omitnan");
    datastd = std(stimeseriesdata,1,1,'omitnan');
end

% datastd = stdmean(datamean,)
% set up figure
if isempty(hfig) || ~ishandle(hfig),
  hfig = figure('Visible',visible);
end

ts = trx.movie_timestamps{1}(1:nframes);
clf(hfig);
set(hfig,'Units','pixels','Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfig);

% plot the indicator patches
f0s = ind.starton;
f1s = ind.endon;

t0s = trx.movie_timestamps{1}(f0s);
t1s = trx.movie_timestamps{1}(f1s);
x_ind = [t0s;t1s;t1s;t0s];
y_lim(1) = min(min(datamean-datastd),-.5);
y_lim(2) = max(datamean+datastd,[],'omitnan');
y_ind = [y_lim(2);y_lim(2);y_lim(1);y_lim(1)];
y_ind = repmat(y_ind,1,numel(t0s));

%plot indicators in back 
patch(hax,x_ind,y_ind,colors.stimbkgd,'LineStyle','none')
hold(hax,"on")

% assume only ALL nan data at the end of movie
ndata = sum(~isnan(datamean));
off0s = [1,f1s];
off1s = [f0s,ndata];
%check
assert(sum(isnan(datamean(1:ndata))) == 0)

if ~ispercbeh
       %plot the std patches during stims
    for i = 1:numel(f0s)
        x = ts(f0s(i):f1s(i));
        yneg = datamean(f0s(i):f1s(i))-datastd(f0s(i):f1s(i));
        ypos = datamean(f0s(i):f1s(i))+datastd(f0s(i):f1s(i));
        patch(hax,[x, fliplr(x)],[yneg,fliplr(ypos)],colors.stimstd,'LineStyle','none')
        hold(hax,"on")
    end
    %plot the std patches during stims
    for i = 1:numel(off0s)
        x = ts(off0s(i):off1s(i));
        yneg = datamean(off0s(i):off1s(i))-datastd(off0s(i):off1s(i));
        ypos = datamean(off0s(i):off1s(i))+datastd(off0s(i):off1s(i));
        patch(hax,[x, fliplr(x)],[yneg,fliplr(ypos)],colors.offstd,'LineStyle','none')
        hold(hax,"on")
    end
    % plot during stim means
    for i = 1:numel(f0s)
        x = ts(f0s(i):f1s(i));
        y = datamean(f0s(i):f1s(i));
        plot(hax,x,y,'-','LineWidth',1,'Color',colors.stim);
        hold(hax,"on")
    end
    %plot during not stim means
    for i = 1:numel(off0s)
        x = ts(off0s(i):off1s(i));
        y = datamean(off0s(i):off1s(i));
        plot(hax,x,y,'-','LineWidth',1,'Color',colors.off);
        hold(hax,"on")
    end
else

    for i = 1:numel(f0s)
        x = ts(f0s(i):f1s(i));
        y = datamean(f0s(i):f1s(i));
        plot(hax,x,y,'-','LineWidth',1,'Color',colors.stim);
    end
    %plot during not stim means
    for i = 1:numel(off0s)
        x = ts(off0s(i):off1s(i));
        y = datamean(off0s(i):off1s(i));
        plot(hax,x,y,'-','LineWidth',1,'Color',colors.off);
    end
end

%set figure labels

if ispercbeh
    [~,behname] = regexp(field,'scores_','match','split');
    ylabel(hax,sprintf('%% of flies performing %s',behname{2}),'Interpreter','none','Fontsize',fontsize);
    title(hax,sprintf('Time series of %s',behname{2}),'Interpreter','none');
    y_lim(1)= -5;
else
    ylabel(hax,field,'Interpreter','none','Fontsize',fontsize);
    title(hax,sprintf('Time series of %s',field),'Interpreter','none');
    y_lim(1) = min(datamean-datastd);
end
xlabel(hax,'Time (s)','Fontsize',fontsize);
set(hax,'YLim',y_lim);

