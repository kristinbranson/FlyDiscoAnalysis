function [phasediff_hilbert] = computeContinuousPhaseDiff_hilbert(norm_ytips,loctall,locball, currwalk_tips_pos_body_Y,w,debug)

nlimb = size(norm_ytips,1);
nwalkfrms = size(norm_ytips,2);

phasediff_hilbert = struct;
legphases = nan(nlimb,nwalkfrms);


for limb = 1:nlimb
    loct = loctall{limb}(2,:);
    locb = locball{limb}(2,:);
    % compute per leg phases with the hilbert transform
    data = angle(hilbert(norm_ytips(limb,:)));

    sidx = min(loct(1),locb(1));
    eidx = max(loct(end),locb(end));
    % only use data after first peak/trough and before last peak/trough
    legphases(limb,sidx:eidx) = data(sidx:eidx);

end

if debug
    % test plots
    % NOTE using diff so limbs(2) - limbs(1);
    % matched code below - make first entry reference 
    limbs = [3,5];
    if ~isempty(legphases(~isnan(legphases(limbs(1),:)))) & ~isempty(legphases(~isnan(legphases(limbs(2),:))))

        figure('Position',[400 747 827 883])
        tiledlayout(6, 1);

        nexttile;
        plot(currwalk_tips_pos_body_Y(limbs,:)');
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('bodyref Y of tips');
        title(sprintf('walk: %d local hilbert phase analysis of limb %d minus limb %d',w,limbs(1),limbs(2)))
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
        limbref = [limbs(2),limbs(1)];
        plot(wrapTo2Pi(diff(unwrap(legphases(limbref,:),[],2))'),'.')
        set(gca,'Ylim',[0,2*pi]);
        set(gca,'XLim',[0,size(norm_ytips,2)+5])
        ylabel('phase diff')

        nexttile
        polarhistogram(diff(legphases(limbref,:))')

    else
        sprintf('walk no. %d: \nlimb %d had %d top peaks and %d bottom peaks\nlimb %d had %d top peaks and %d bottom peaks', ...
            w, limbs(1), size(loctall{limbs(1)},2), size(locball{limbs(1)},2),...
            limbs(2), size(loctall{limbs(2)},2), size(locball{limbs(2)},2))
    end
end


% compute diffs
legorder = {'RF','RM','RH','LH','LM','LF'};

% first leg is reference
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
    2,3
    1,5
    6,2
    2,4];

for d = 1:numel(diffs)/2
    name = [legorder{diffs(d,1)},'_',legorder{diffs(d,2)}];
    % NOTE using diff so limbs(2) - limbs(1) so reversing diffs values to
    % make reference correct.
    % limbs = [diffs(d,2),diffs(d,1)]; % limbs was reveresed order when using diff
    % currdata1 = wrapTo2Pi(diff(unwrap(legphases(limbs,:),[],2)));
        % currdata1 = wrapToPi(diff(unwrap(legphases(limbs,:),[],2)));

    % compute circ_dist instead ... 
    limbs = [diffs(d,1),diffs(d,2)];
    currdata = circ_dist(legphases(limbs(1),:),legphases(limbs(2),:));

% debug
figure, 
hold on
plot(legphases(limbs(1),:),'.k')
plot(legphases(limbs(2),:),'.r')

    phasediff_hilbert.(name).data = currdata;
    % reported in -pi to pi
    phasediff_hilbert.(name).mean = circ_mean(currdata(~isnan(currdata)));
    phasediff_hilbert.(name).std = circ_std(currdata(~isnan(currdata)));
    phasediff_hilbert.(name).n = numel(currdata);

end

% group ipsi and contra
phasegroups = struct;
phasegroups.ipsi_post_2 = {'LM_LH','RM_RH'};
phasegroups.ipsi_ant_2 = {'LF_LM','RF_RM'};
phasegroups.ipsi_A2P_4 = {'LM_LH','RM_RH','LF_LM','RF_RM'};
phasegroups.conta_L2R_3 = {'LF_RF','LM_RM','LH_RH'};
phasegroups.tripods_4 = {'RF_LM','LM_RH','LF_RM','RM_LH'};


% combine phase groups
flds = fields(phasegroups);
for f= 1:numel(flds)
phasegroupname = flds{f};
phasegroup = phasegroups.(phasegroupname);

curr_all = {};
for p = 1:numel(phasegroup)
    curr_all{p} = phasediff_hilbert.(phasegroup{p}).data;
end
diffdata = [];
diffdata = horzcat(curr_all{:});
phasediff_hilbert.(phasegroupname).data = diffdata;
phasediff_hilbert.(phasegroupname).mean = circ_mean(diffdata(~isnan(diffdata)));
phasediff_hilbert.(phasegroupname).std = circ_std(diffdata(~isnan(diffdata)));
phasediff_hilbert.(phasegroupname).n = numel(diffdata(~isnan(diffdata)));


end
% if debug
% figure
% hold on
% tiledlayout(3, 2);
% 
% for iii = 1:numel(flds)
% nexttile
% polarhistogram(phasediff_hilbert.(flds{iii}).diffdata)
% end 
% 
% 
% end


phasediff_hilbert.phasedata = legphases;

end