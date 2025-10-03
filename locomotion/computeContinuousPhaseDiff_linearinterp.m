function [phasediff_interp] = computeContinuousPhaseDiff_linearinterp(norm_ytips,currwalk_tips_pos_body_Y,w)

debug = 1;
minPeakProminence = 1;
nlimb = size(norm_ytips,1);
nwalkfrms = size(norm_ytips,2);

phasediff_interp = struct;
legphases = nan(nlimb,nwalkfrms);
loctall = {};
locball = {};
for limb = 1:nlimb
    [toppeaks, loct] = findpeaks(norm_ytips(limb,:),'MinPeakProminence', minPeakProminence);
    [bottompeaks, locb] = findpeaks((1-norm_ytips(limb,:)),'MinPeakProminence', minPeakProminence);
    loctall{limb}  = [toppeaks;loct];
    locball{limb} = [bottompeaks; locb];

    % step defined as swing (min to max) + stance (max to min) but
    % don't want to loose half steps at beginning or end





    if numel(locb) >= 3 & numel(loct) >= 3
        % if walk starts with swing for this leg
        % pair trough-peak-trough, skip if no pair, two peaks
        for l = 1:numel(locb)-1
            idx_peaks = find(loct < locb(l+1) & loct > locb(l));
            peak_loc = loct(idx_peaks);
            if numel(peak_loc) == 1
                tipdata_swing = norm_ytips(limb,locb(l):peak_loc);
                tipdata_stance = norm_ytips(limb,peak_loc:locb(l+1));
                swingphase = pi * (tipdata_swing-tipdata_swing(1)) / (tipdata_swing(end) - tipdata_swing(1));
                stancephase =  (pi * (tipdata_stance-tipdata_stance(1)) / (tipdata_stance(end) - tipdata_stance(1)))+ pi;
                legphases(limb,locb(l):locb(l+1)) = [swingphase,stancephase(2:end)];


            end

        end
        % if walk starts with stance phase for this leg
        if loct(1) < locb(1) & numel(loct) >= 3
            % add the first stance to phases
            tipdata_stance = norm_ytips(limb,loct(1):locb(1));
            stancephase =  (pi * (tipdata_stance-tipdata_stance(1)) / (tipdata_stance(end) - tipdata_stance(1)))+ pi;
            legphases(limb,loct(1):locb(1)) = stancephase;
        end
        % if walk ends with swing for this leg
        if loct(end) > locb(end)
            tipdata_swing = norm_ytips(limb,locb(end):loct(end));
            swingphase = pi * (tipdata_swing-tipdata_swing(1)) / (tipdata_swing(end) - tipdata_swing(1));
            legphases(limb,locb(end):loct(end)) = swingphase;


        end
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
        title(sprintf('walk %d interp phase analysis of limb %d vs limb %d',w, limbs(2),limbs(1)))
        legend(sprintf('limb %d',limbs(1)),sprintf('limb %d',limbs(2)))
        nexttile;
        hold on
        % mean centered and normalized data
        plot(norm_ytips(limbs,:)');
        % plot detected peaks
        plot(loctall{limbs(1)}(2,:),loctall{limbs(1)}(1,:),'ob')
        plot(locball{limbs(1)}(2,:),1-locball{limbs(1)}(1,:),'xb')
        plot(loctall{limbs(2)}(2,:),loctall{limbs(2)}(1,:),'or')
        plot(locball{limbs(2)}(2,:),1-locball{limbs(2)}(1,:),'xr')
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('zscore y pos, pks');
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
            w, limbs(1), size(loctall{limbs(1)},2), size(locball{limbs(1)},2),...
            limbs(2), size(loctall{limbs(2)},2), size(locball{limbs(2)},2))

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