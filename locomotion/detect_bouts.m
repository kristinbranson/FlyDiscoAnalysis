%AR 20240602 converted from python code with chatGPT assist
function [start_indices,end_indices] = detect_bouts(data, include_edge_bouts)
    if nargin < 2
        include_edge_bouts = false;
    end
    assert(isvector(data), 'Input data must be a 1D array');
    data = double(data(:));
    transitions = diff(data);

    % starts are transitions from 0 to 1, add 1 to account for diff
        % offset
        start_indices = find(transitions == 1) + 1;
         % ends are transitions from 1 to 0
        end_indices = find(transitions == -1) + 1;
        if include_edge_bouts
            % Include bouts ongoing at signal start/end
            if data(1) == 1
                start_indices = [1; start_indices];
            end
            if data(end) == 1
                end_indices = [end_indices; numel(data) + 1];
            end
        else
            % Default: discard incomplete bouts at signal edges
            if data(1) == 1
                end_indices = end_indices(2:end);
            end
            if data(end) == 1
                start_indices = start_indices(1:end-1);
            end
        end

end