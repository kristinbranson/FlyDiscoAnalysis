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
pxpermm = trx.pxpermm;
bout_metrics = struct;
% list of perframe features to compute stats over bouts
pfflist = {'velmag_ctr','du_ctr','dv_ctr','dtheta','absdtheta','CoM_stability'};
% perframe features if 'first' = 1st derivative (d) 'second' = second
% derivative (dd), or 'none' = no derivative
pfflist_flags = {'first','first','first','first','first','none'};


%compute metrics for each limb of each fly

for is = 1:numel(state)
    %     is
    for fly = 1:numel(boutstruct)
        %         fly
        currfly_timestamps = timestamps(trx.firstframes(fly):trx.endframes(fly));
        for limb = 1:numel(boutstruct(fly).perlimb)
            %             limb
            if ismember(state{is},{'swing','stance'})

                start_indices = boutstruct(fly).perlimb(limb).(state{is}).start_indices;
                end_indices = boutstruct(fly).perlimb(limb).(state{is}).end_indices;


                % compute frame and time duration (ms) from timestamps
                [durations_frames,durations_time] = computeBoutDurations(start_indices,end_indices,currfly_timestamps);
                %durations
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).durations_frames = durations_frames;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).durations_time = durations_time;

                % % made more general AR 20250918

                % loop over list of perframe features
                for ifns = 1:numel(pfflist)
                    fn = pfflist{ifns};
                    derivative_flag = pfflist_flags{ifns};
                    datastruct = compute_StatsofPreframeFeatureDuringBouts(fly,fn,trx,start_indices,end_indices,derivative_flag);
                    bout_metrics.perfly(fly).perlimb(limb).(state{is}).(fn) = datastruct;
                end

                % compute leg tip speeds in body coordinates AR 9/18/25 TO DO make general for legperframefeatures
                % speed in mm/s
                limb_idx = limb;
                [mean_tips_speed_bodyref,std_tips_speed_bodyref,min_speed_bodyref,max_speed_bodyref] = compute_legtipspeed(fly,tips_pos_body,limb_idx,start_indices,end_indices,currfly_timestamps,pxpermm);
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).mean_tips_speed_bodyref = mean_tips_speed_bodyref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).std_tips_speed_bodyref = std_tips_speed_bodyref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).min_speed_bodyref = min_speed_bodyref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).max_speed_bodyref = max_speed_bodyref;

                % compute leg tip speeds in global coordinates
                position_data = aptdata.pTrk;
                limb_idx = legtip_landmarknums(limb);
                [mean_tips_speed_globalref,std_tips_speed_globalref,min_speed_globalref,max_speed_globalref] = compute_legtipspeed(fly,position_data,limb_idx,start_indices,end_indices,currfly_timestamps,pxpermm);
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).mean_tips_speed_globalref = mean_tips_speed_globalref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).std_tips_speed_globalref = std_tips_speed_globalref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).min_speed_globalref = min_speed_globalref;
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).max_speed_globalref = max_speed_globalref;

                % pass along start and end indices
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).start_indices = start_indices';
                bout_metrics.perfly(fly).perlimb(limb).(state{is}).end_indices = end_indices';

            end

            if strcmp(state{is},'step')

                step_t0s = boutstruct(fly).perlimb(limb).(state{is}).start_indices;
                step_t1s = boutstruct(fly).perlimb(limb).(state{is}).end_indices;
                stance_t0s = boutstruct(fly).perlimb(limb).stance.start_indices;
                stance_t1s = boutstruct(fly).perlimb(limb).stance.end_indices;
                tip_pos_body = tips_pos_body{fly};


                stepfeatures = computeStepFeatures(fly,trx,aptdata,tip_pos_body,legtip_landmarknums,limb,step_t0s,step_t1s,stance_t0s,stance_t1s,currfly_timestamps);

                bout_metrics.perfly(fly).perlimb(limb).(state{is}) = stepfeatures;

                % compute mean and standard deviation for each boutfeatures


            end



        end
        % combine into pairs and all_limbs for each fly
        flds = fields(bout_metrics.perfly(fly).perlimb(limb).(state{is}));

        for fld = 1:numel(flds)
            % pairs
            for pair = 1:size(pairs,1)
                if ~isstruct(bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}))
                    bout_metrics.perfly(fly).pairs(pair).(state{is}).(flds{fld}) = [bout_metrics.perfly(fly).perlimb(pairs(pair,1)).(state{is}).(flds{fld}), ...
                        bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld})];
                elseif isstruct(bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}))
                    subflds = fields(bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}));
                    for sbf = 1:numel(subflds)

                        bout_metrics.perfly(fly).pairs(pair).(state{is}).(flds{fld}).(subflds{sbf}) = [bout_metrics.perfly(fly).perlimb(pairs(pair,1)).(state{is}).(flds{fld}).(subflds{sbf}), ...
                            bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}).(subflds{sbf})];
                    end
                end
            end
            % all_limbs
            if ~isstruct(bout_metrics.perfly(fly).perlimb(1).(state{is}).(flds{fld}))
                currfly_all_limbs = {};
                for limb = 1:numel(boutstruct(fly).perlimb)
                    currfly_all_limbs{limb} = bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld});
                end
                bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}) = horzcat(currfly_all_limbs{:});
            % perframe feature structs loop over the fields of the
            % datastruct
            elseif isstruct(bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld}))

                subflds = fields(bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}));
                for sbf = 1:numel(subflds)
                    currfly_all_limbs = {};
                    for limb = 1:numel(boutstruct(fly).perlimb)
                        currfly_all_limbs{limb} = bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld}).(subflds{sbf});
                    end
                    bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}).(subflds{sbf}) = horzcat(currfly_all_limbs{:});
                end
            end
        end
    end
end



% combine perfly to allflies for perlimb, pairs, and all_limbs
% TO DO - more sophfisticated way to combine data across flies (extend to
% across movies?)

for is = 1:numel(state)
    flds = fields(bout_metrics.perfly(1).perlimb(1).(state{is}));
    for fld = 1:numel(flds)
        
        if ~isstruct(bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}))
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
        elseif isstruct(bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}))
            subflds = fields(bout_metrics.perfly(fly).perlimb(pairs(pair,2)).(state{is}).(flds{fld}));
            for sbf = 1:numel(subflds)
                % all limbs
                curr_allfly=  {};
                for fly = 1:numel(bout_metrics.perfly)
                    curr_allfly{fly} = (bout_metrics.perfly(fly).all_limbs.(state{is}).(flds{fld}).(subflds{sbf}));
                end
                bout_metrics.allflies.all_limbs.(state{is}).(flds{fld}).(subflds{sbf}) = horzcat(curr_allfly{:});
                % per limb
                for limb = 1:nlimb
                    curr_allfly = {};
                    for fly = 1:nflies
                        curr_allfly{fly} = (bout_metrics.perfly(fly).perlimb(limb).(state{is}).(flds{fld}).(subflds{sbf}));
                    end
                    bout_metrics.allflies.perlimb(limb).(state{is}).(flds{fld}).(subflds{sbf}) = horzcat(curr_allfly{:});
                end
                % per pairs

                for pair = 1:npairs
                    curr_allfly = {};
                    for fly = 1:nflies
                        curr_allfly{fly} = (bout_metrics.perfly(fly).pairs(pair).(state{is}).(flds{fld}).(subflds{sbf}));
                    end
                    bout_metrics.allflies.pairs(pair).(state{is}).(flds{fld}).(subflds{sbf}) = horzcat(curr_allfly{:});
                end
            end
        end
    end
end

% updated to match general schema change above - AR 20250918
% compute mean bout durations in speed bins
% % only run for swing and stance
for is = 1:2
    for pair = 1:npairs
        [Nboutspervelmagbin,velmagbincenters,meanboutdurationsofvelmagbins,stdboutdurationsofvelmagbins] = compute_meanboutdurationVSspeed(bout_metrics.allflies.pairs(pair).(state{is}).durations_time,...
            bout_metrics.allflies.pairs(pair).(state{is}).velmag_ctr.mean,'binedges',binedges);
        bout_metrics.allflies.pairs(pair).(state{is}).Nboutspervelmagbin =Nboutspervelmagbin;
        bout_metrics.allflies.pairs(pair).(state{is}).velmagbincenters = velmagbincenters;
        bout_metrics.allflies.pairs(pair).(state{is}).meanboutdurationsofvelmagbins = meanboutdurationsofvelmagbins;
        bout_metrics.allflies.pairs(pair).(state{is}).stdboutdurationsofvelmagbins = stdboutdurationsofvelmagbins;
    end
    % for all limbs

    [Nboutspervelmagbin,velmagbincenters,meanboutdurationsofvelmagbins,stdboutdurationsofvelmagbins] = compute_meanboutdurationVSspeed(bout_metrics.allflies.all_limbs.(state{is}).durations_time,...
        bout_metrics.allflies.all_limbs.(state{is}).velmag_ctr.mean,'binedges',binedges);
    bout_metrics.allflies.all_limbs.(state{is}).Nboutspervelmagbin =Nboutspervelmagbin;
    bout_metrics.allflies.all_limbs.(state{is}).velmagbincenters = velmagbincenters;
    bout_metrics.allflies.all_limbs.(state{is}).meanboutdurationsofvelmagbins = meanboutdurationsofvelmagbins;
    bout_metrics.allflies.all_limbs.(state{is}).meanboutdurationsofvelmagbins = stdboutdurationsofvelmagbins;

end



