function [bout_metrics] = computeboutmetrics2(trx,aptdata,tips_pos_body,legtip_landmarknums,boutstruct)

% input struct is output of detect_bouts or find_bout_overlap
% name(fly).perlimb(1:6).start_indices
% name(fly).perlimb(1:6).end_indices

%output fields 


state = {'swing','stance','step'};
pairs = [1,6;2,5;3,4];
binedges = 5:2:34;
pairNames = {'Front','Mid','Rear'};
npairs = size(pairs,1);
nlimb = numel(boutstruct(1).perlimb);
nflies = numel(boutstruct);
timestamps = trx.movie_timestamps{1};
bout_metrics = struct;


%compute metrics for each limb of each fly

for is = 1:numel(state)
%     is
    for fly = 1:numel(boutstruct)
%         fly
        timestamps_perfly = timestamps(trx.firstframes(fly):trx.endframes(fly));
        for limb = 1:numel(boutstruct(fly).perlimb)
%             limb
            if ismember(state{is},{'swing','stance'})
            % compute frame and time duration
            
            [durations_frames,durations_time] = computeBoutDurations(boutstruct(fly).perlimb(limb).(state{is}).start_indices, ...
                boutstruct(fly).perlimb(limb).(state{is}).end_indices,timestamps_perfly);

            % compute mean body speed during bout
            [meanvelmag,stdnvelmag] = computeMeanPreframeDuringBout(fly,'velmag',trx,boutstruct(fly).perlimb(limb).(state{is}).start_indices,boutstruct(fly).perlimb(limb).(state{is}).end_indices);

            %compile data
            %durations
            bout_metrics.perfly(fly).perlimb(limb).(state{is}).durations_frames = durations_frames;
            bout_metrics.perfly(fly).perlimb(limb).(state{is}).durations_time = durations_time;


            %perframe values
            bout_metrics.perfly(fly).perlimb(limb).(state{is}).meanvelmag = meanvelmag;
            bout_metrics.perfly(fly).perlimb(limb).(state{is}).stdvelmag = stdnvelmag;

            end
            if strcmp(state{is},'step')

            step_t0s = boutstruct(fly).perlimb(limb).(state{is}).start_indices;
            step_t1s = boutstruct(fly).perlimb(limb).(state{is}).end_indices;
            stance_t0s = boutstruct(fly).perlimb(limb).stance.start_indices;
            stance_t1s = boutstruct(fly).perlimb(limb).stance.end_indices;
            tip_pos_body = tips_pos_body{fly};
            

            boutfeatures = computeBoutFeatures(fly,trx,aptdata,tip_pos_body,legtip_landmarknums,limb,step_t0s,step_t1s,stance_t0s,stance_t1s,timestamps_perfly);

            bout_metrics.perfly(fly).perlimb(limb).(state{is}) = boutfeatures;

            % compute mean and standard deviation for each boutfeatures

               
            end



        end
        % combine into pairs and all_limbs for each fly
        flds = fields(bout_metrics.perfly(fly).perlimb(limb).(state{is}));

        for fld = 1:numel(flds)
            % pairs
            for pair = 1:size(pairs,1)
                bout_metrics.perfly(fly).pairs(pair).(state{is}).(flds{fld}) = [bout_metrics.perfly(fly).perlimb(pairs(pair,1)).(state{is}).(flds{fld}), ...
                    bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld})];
            end
            % all_limbs
            currfly_all_limbs = {};
            for limb = 1:numel(boutstruct(fly).perlimb)
                currfly_all_limbs{limb} = bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld});
            end
            bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}) = horzcat(currfly_all_limbs{:});

        end
    end
end



% combine perfly to allflies for perlimb, pairs, and all_limbs


for is = 1:numel(state)
    flds = fields(bout_metrics.perfly(1).perlimb(1).(state{is}));
    for fld = 1:numel(flds)

        % all limbs
        curr_allfly=  {};
        for fly = 1:numel(bout_metrics.perfly)
            curr_allfly{fly} = (bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}));
        end
        bout_metrics.allflies.all_limbs.(state{is}).(flds{fld}) = horzcat(curr_allfly{:});

        % per limb

        for limb = 1:nlimb
            curr_allfly = {};
            for fly = 1:nflies
                curr_allfly{fly} = (bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld}));
            end
            bout_metrics.allflies.perlimb(limb).(state{is}).(flds{fld}) = horzcat(curr_allfly{:});
        end

        % per pairs

        for pair = 1:npairs
            curr_allfly = {};
            for fly = 1:nflies
                curr_allfly{fly} = (bout_metrics.perfly(fly).pairs(pair).(state{is}).(flds{fld}));
            end
            bout_metrics.allflies.pairs(pair).(state{is}).(flds{fld}) = horzcat(curr_allfly{:});
        end
    end
end

% compute mean bout durations in speed bins
% only run for swing and stance 
for is = 1:2
    for pair = 1:npairs
        [Nboutspervelmagbin,velmagbincenters,meanboutdurationsofvelmagbins,stdboutdurationsofvelmagbins] = compute_meanboutdurationVSspeed(bout_metrics.allflies.pairs(pair).(state{is}).durations_time,...
            bout_metrics.allflies.pairs(pair).(state{is}).meanvelmag,'binedges',binedges);
        bout_metrics.allflies.pairs(pair).(state{is}).Nboutspervelmagbin =Nboutspervelmagbin;
        bout_metrics.allflies.pairs(pair).(state{is}).velmagbincenters = velmagbincenters;
        bout_metrics.allflies.pairs(pair).(state{is}).meanboutdurationsofvelmagbins = meanboutdurationsofvelmagbins;
        bout_metrics.allflies.pairs(pair).(state{is}).stdboutdurationsofvelmagbins = stdboutdurationsofvelmagbins;
    end
    % for all limbs
    
    [Nboutspervelmagbin,velmagbincenters,meanboutdurationsofvelmagbins,stdboutdurationsofvelmagbins] = compute_meanboutdurationVSspeed(bout_metrics.allflies.all_limbs.(state{is}).durations_time,...
        bout_metrics.allflies.all_limbs.(state{is}).meanvelmag,'binedges',binedges);
    bout_metrics.allflies.all_limbs.(state{is}).Nboutspervelmagbin =Nboutspervelmagbin;
    bout_metrics.allflies.all_limbs.(state{is}).velmagbincenters = velmagbincenters;
    bout_metrics.allflies.all_limbs.(state{is}).meanboutdurationsofvelmagbins = meanboutdurationsofvelmagbins;
    bout_metrics.allflies.all_limbs.(state{is}).meanboutdurationsofvelmagbins = stdboutdurationsofvelmagbins;

end



