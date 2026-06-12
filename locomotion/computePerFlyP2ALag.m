function perflylag = computePerFlyP2ALag(walkfeaturestruct, fieldname)
% computePerFlyP2ALag - Aggregate per-walk P2A onset lag into per-fly.
%
% Generic over the metric field (e.g. 'Pliftoff2Aliftoff_lag' or
% 'Ptouchdown2Aliftoff_lag'). Each per-walk entry has fieldname with the 7
% sub-fields (RH_to_RM, RM_to_RF, LH_to_LM, LM_to_LF, H_to_M, M_to_F, all),
% each with .data (per-event lags, ms) and .mean.
%
% Dual aggregation (matching computePerFlyphasefeatures, but linear stats):
%   frm_* - concatenate per-event .data across walks, then stats
%   walk_* - per-walk .mean values, then stats
%
% Output: perflylag.(subfield) with
%   .frm_data_fly/.frm_mean_fly/.frm_std_fly/.frm_n_fly
%   .walk_data_fly/.walk_mean_fly/.walk_std_fly/.walk_n_fly

nwalks = numel(walkfeaturestruct);
subfields = fieldnames(walkfeaturestruct(1).(fieldname));

for sf = 1:numel(subfields)
    sn = subfields{sf};

    frm_cells = cell(1, nwalks);
    walk_means = nan(1, nwalks);
    for w = 1:nwalks
        frm_cells{w} = walkfeaturestruct(w).(fieldname).(sn).data;
        walk_means(w) = walkfeaturestruct(w).(fieldname).(sn).mean;
    end

    frmdata = horzcat(frm_cells{:});
    perflylag.(sn).frm_data_fly = frmdata;
    perflylag.(sn).frm_mean_fly = mean(frmdata, 'omitnan');
    perflylag.(sn).frm_std_fly  = std(frmdata, 'omitnan');
    perflylag.(sn).frm_n_fly    = sum(~isnan(frmdata));

    perflylag.(sn).walk_data_fly = walk_means;
    perflylag.(sn).walk_mean_fly = mean(walk_means, 'omitnan');
    perflylag.(sn).walk_std_fly  = std(walk_means, 'omitnan');
    perflylag.(sn).walk_n_fly    = sum(~isnan(walk_means));
end

end
