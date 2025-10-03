function [phasediff_interp] = computeContinuousPhaseDiff_hilbert(norm_ytips,currwalk_tips_pos_body_Y,w)

debug = 1;
minPeakProminence = 1;
nlimb = size(norm_ytips,1);
nwalkfrms = size(norm_ytips,2);

phasediff_interp = struct;
legphases = nan(nlimb,nwalkfrms);


for limb = 1:nlimb
    [~, loct] = findpeaks(norm_ytips(limb,:),'MinPeakProminence', minPeakProminence);
    [~, locb] = findpeaks((1-norm_ytips(limb,:)),'MinPeakProminence', minPeakProminence);
        loctall{limb}  = [loct];
    locball{limb} = [locb];
    data = angle(hilbert(norm_ytips(limb,:)));
    if numel(locb) >= 3 & numel(loct) >= 3
        sidx = min(loct(1),locb(1));
        eidx = max(loct(end),locb(end));
        % only use data after first peak/trough and before last peak/trough
        legphases(limb,sidx:eidx) = data(sidx:eidx);
    end
end

if debug
    % test plots
    % NOTE using diff so limbs(2) - limbs(1);
    limbs = [1,3];
    if ~isempty(legphases(~isnan(legphases(limbs(1),:)))) & ~isempty(legphases(~isnan(legphases(limbs(2),:))))

        figure('Position',[1326 747 827 883])

        tiledlayout(6, 1);
        nexttile;
        plot(currwalk_tips_pos_body_Y(limbs,:)');
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('bodyref Y of tips');
        title(sprintf('walk: %d hilbert phase analysis of limb %d vs limb %d',w,limbs(2),limbs(1)))
        legend(sprintf('limb %d',limbs(1)),sprintf('limb %d',limbs(2)))

        nexttile;
        hold on
        % mean centered and normalized data
        plot(norm_ytips(limbs,:)');
        ylabel('zscore y pos');
        set(gca,'XLim',[0,size(norm_ytips,2)+5])

        nexttile([2,1]);
        plot(1:size(norm_ytips,2),(legphases(limbs,:)'),'.')
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('phases')

        nexttile
        plot(wrapTo2Pi(diff(unwrap(legphases(limbs,:),[],2))'),'.')
        set(gca,'Ylim',[0,2*pi]);
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('phase diff')

        nexttile
        polarhistogram(diff(legphases(limbs,:))')

    else
        sprintf('walk no. %d: \nlimb %d had %d top peaks and %d bottom peaks\nlimb %d had %d top peaks and %d bottom peaks', ...
            w, limbs(1), numel(loctall{limbs(1)}), numel(locball{limbs(1)}),...
            limbs(2), numel(loctall{limbs(2)}), numel(locball{limbs(2)}))

    end
end


% compute diffs
legorder = {'RF','RM','RH','LH','LM','LF'};

% first leg in reference
diffs = [5, 1
    5, 2
    5, 3
    5, 4
    5, 6
    6,1
    5,2
    4,3
    6,5
    1,2
    2,3];

for d = 1:numel(diffs)/2
    name = [legorder{diffs(d,1)},'_',legorder{diffs(d,2)}];
    % NOTE using diff so limbs(2) - limbs(1) so reversing diffs values to
    % make reference correct.
    limbs = [diffs(d,2),diffs(d,1)];
    currdata = wrapTo2Pi(diff(unwrap(legphases(limbs,:),[],2)));
    phasediff_interp.(name).data = currdata;
    % reported in -pi to pi
    phasediff_interp.(name).mean = circ_mean(currdata(~isnan(currdata)));
    phasediff_interp.(name).std = circ_std(currdata(~isnan(currdata)));


end
end