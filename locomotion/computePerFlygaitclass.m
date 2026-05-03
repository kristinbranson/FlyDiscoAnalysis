function perflygait = computePerFlygaitclass(walkfeaturestruct)
% computePerFlygaitclass - Aggregate per-walk gait_class structs into per-fly.
%
% Input:
%   walkfeaturestruct - struct array, each element has .gait_class with
%                       .data (1 x nframes class codes) and per-class
%                       counts (tripod_count, tetrapod_count, ...).
%
% Output:
%   perflygait - struct with:
%     .frm_data_fly        - concatenated class codes across all walks (uint8)
%     .frm_n_fly           - total frame count
%     .tripod_count_fly    - sum across walks
%     .tetrapod_count_fly
%     .grounded_count_fly
%     .airborne_count_fly
%     .other_count_fly
%
% Class codes per gait_pattern_constants.m: 1=tripod 2=tetrapod
% 3=grounded 4=airborne 5=other.

class_names = {'tripod', 'tetrapod', 'grounded', 'airborne', 'other'};

nwalks = numel(walkfeaturestruct);

% Concatenate raw per-frame class codes across walks
data_cells = cell(1, nwalks);
for w = 1:nwalks
    data_cells{w} = walkfeaturestruct(w).gait_class.data;
end
perflygait.frm_data_fly = horzcat(data_cells{:});
perflygait.frm_n_fly    = numel(perflygait.frm_data_fly);

% Sum per-class counts across walks
for c = 1:numel(class_names)
    cf = [class_names{c} '_count'];
    counts = zeros(1, nwalks);
    for w = 1:nwalks
        counts(w) = walkfeaturestruct(w).gait_class.(cf);
    end
    perflygait.([class_names{c} '_count_fly']) = sum(counts);
end

end
