function perexplag = computePerExpP2ALag(perwalk_metrics, perfly_metrics, fieldname)
% computePerExpP2ALag - Aggregate per-walk P2A onset lag into per-experiment.
%
% Event-pooled across all walks (equivalently all flies) in the experiment,
% mirroring computePerExpTCS. Generic over the metric field (e.g.
% 'Pliftoff2Aliftoff_lag' or 'Ptouchdown2Aliftoff_lag').
%
% Inputs:
%   perwalk_metrics - struct array (one entry per walk across all flies),
%                     each with fieldname and its 7 sub-fields.
%   perfly_metrics  - struct array (one entry per fly); unused, kept to
%                     mirror computePerExpTCS' signature.
%   fieldname       - metric field name.
%
% Output: perexplag.(subfield) with
%   .frm_data_exp/.frm_mean_exp/.frm_std_exp/.frm_n_exp
%   .walk_data_exp/.walk_mean_exp/.walk_std_exp/.walk_n_exp

nwalks = numel(perwalk_metrics);
subfields = fieldnames(perwalk_metrics(1).(fieldname));

for sf = 1:numel(subfields)
    sn = subfields{sf};

    frm_cells = cell(1, nwalks);
    walk_means = nan(1, nwalks);
    for w = 1:nwalks
        frm_cells{w} = perwalk_metrics(w).(fieldname).(sn).data;
        walk_means(w) = perwalk_metrics(w).(fieldname).(sn).mean;
    end

    frmdata = horzcat(frm_cells{:});
    perexplag.(sn).frm_data_exp = frmdata;
    perexplag.(sn).frm_mean_exp = mean(frmdata, 'omitnan');
    perexplag.(sn).frm_std_exp  = std(frmdata, 'omitnan');
    perexplag.(sn).frm_n_exp    = sum(~isnan(frmdata));

    perexplag.(sn).walk_data_exp = walk_means;
    perexplag.(sn).walk_mean_exp = mean(walk_means, 'omitnan');
    perexplag.(sn).walk_std_exp  = std(walk_means, 'omitnan');
    perexplag.(sn).walk_n_exp    = sum(~isnan(walk_means));
end

end
