function [N,bincenters,meanboutdurationsATspeed,stdboutdurationsATspeed] = compute_meanboutdurationVSspeed(bouts_durations_time,bouts_meanvelmag,varargin)



% Process the args
[binedges] = ...
  myparse(varargin,...
  'binedges',5:1.5:35);

meanboutdurationsATspeed = nan(1,numel(binedges)-1);
stdboutdurationsATspeed = nan(1,numel(binedges)-1);

[N,~,binidx] = histcounts(bouts_meanvelmag,binedges);
for i = 1:numel(binedges)-1
    bincenters(i) = mean(binedges(i:i+1));
    meanboutdurationsATspeed(i) = mean(bouts_durations_time(binidx == i),'omitnan');
    stdboutdurationsATspeed(i) = std(bouts_durations_time(binidx == i),1,'omitnan');
end
