function phase_data = compute_phase_data_peaks(walkscores, tips_pos_body, varargin)
% ANALYZE_PHASE_DATA Analyzes phase offsets between legs during each walking bout
%
% INPUTS:
%   walkscores - cell array of walk scores for each fly
%   tips_pos_body - cell array of leg tip positions in body coordinates
%
% OPTIONAL PARAMETERS (name-value pairs):
%   'minwalkbout' - minimum walk bout length in frames (default: 15)
%   'pretimewindow' - pre-time window as fraction of 2*pi (default: pi/4)
%   'ref_leg' - reference leg number (default: 6)
%   'minPeakHeight' - minimum peak height for findpeaks (default: 2)
%   'minPeakDistance' - minimum peak distance for findpeaks (default: 5)
%   'minPeakProminence' - minimum peak prominence for findpeaks (default: 5)
%   'debug' - enable debug messages (default: false)
%
% OUTPUT:
%   phase_data - struct array containing phase analysis results for each walk bout
%
% EXAMPLE:
%   phase_data = analyze_phase_data(walkscores, tips_pos_body, ...
%                                  'minwalkbout', 20, 'ref_leg', 5);

% Parse input parameters
p = inputParser;
addRequired(p, 'walkscores');
addRequired(p, 'tips_pos_body');
addParameter(p, 'minwalkbout', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'pretimewindow', pi/4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ref_leg', 6, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minPeakHeight', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'minPeakDistance', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minPeakProminence', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'debug', false, @(x) islogical(x) && isscalar(x));

parse(p, walkscores, tips_pos_body, varargin{:});

% Extract parameters
minwalkbout = p.Results.minwalkbout;
pretimewindow = p.Results.pretimewindow;
ref_leg = p.Results.ref_leg;
minPeakHeight = p.Results.minPeakHeight;
minPeakDistance = p.Results.minPeakDistance;
minPeakProminence = p.Results.minPeakProminence;
debug = p.Results.debug;

% Suppress warning
warning('off', 'signal:findpeaks:largeMinPeakHeight');

% Initialize output structure
phase_data = struct('fly', [], 'walkbout', [], 'walkframes', [], ...
                   'peaksperleg', {}, 'closest_peak', [], 'phase_offsets', []);

ct = 0;
nflies = numel(walkscores);

% Determine number of legs from tips_pos_body structure
if ~isempty(tips_pos_body) && ~isempty(tips_pos_body{1})
    nlegs = size(tips_pos_body{1}, 1);
else
    error('tips_pos_body is empty or invalid');
end

% Main processing loop for each fly
for flyi = 1:nflies
    if debug
        fprintf('Fly:%d    %%%% \n', flyi);
    end
    
    % Detect walking bouts
    [walk_t0s, walk_t1s] = detect_bouts(walkscores{flyi});
    
    % Process each walking bout
    for w = 1:numel(walk_t0s)
        curr_legs = cell(1, nlegs);
        curr_data = struct();
        curr_data.closest_peak = [];
        curr_data.phase_offsets = [];
        
        % Check if walk bout is long enough
        if walk_t1s(w) - walk_t0s(w) > minwalkbout
            ct = ct + 1;
            x = walk_t0s(w):walk_t1s(w);
            y = squeeze(tips_pos_body{flyi}(:, 2, walk_t0s(w):walk_t1s(w)));
            
            % Zero-center y data
            y_zeroed = y - mean(y, 2);
            
            % Store basic walk bout information
            curr_data.fly = flyi;
            curr_data.walkbout = w;
            curr_data.walkframes = x;
            
            % Find peaks for each leg
            for legi = 1:nlegs
                [pks, loc] = findpeaks(y_zeroed(legi, :), ...
                                     'MinPeakHeight', minPeakHeight, ...
                                     'MinPeakDistance', minPeakDistance, ...
                                     'MinPeakProminence', minPeakProminence);
                if ~isempty(loc)
                    curr_legs{legi} = x(loc);
                else
                    if debug
                        fprintf('walk %d, no peaks leg:%d \n', w, legi);
                    end
                end
            end
            
            curr_data.peaksperleg = curr_legs;
            phase_data(ct) = curr_data;
            
        else
            if debug
                fprintf('short walk %d, length %d, start frame %d \n', ...
                       w, walk_t1s(w) - walk_t0s(w), walk_t0s(w));
            end
        end
    end
end

% Phase analysis for each valid walk bout
for w = 1:numel(phase_data)
    curr_data = phase_data(w);
    
    % Check if reference leg has enough peaks
    if numel(curr_data.peaksperleg{ref_leg}) > 2
        mperiod = ceil(mean(diff(curr_data.peaksperleg{ref_leg})));
        
        % Loop over steps (peak pairs)
        for s = 1:numel(curr_data.peaksperleg{ref_leg}) - 1
            refpeak = curr_data.peaksperleg{ref_leg}(s);
            period = curr_data.peaksperleg{ref_leg}(s + 1) - refpeak;
            frame_window = (pretimewindow / (2 * pi)) * period;
            
            % Check if period is reasonable
            if period >= mperiod * 1.5
                if debug
                    fprintf('long period fly %d, frame %d \n', curr_data.fly, refpeak);
                end
            else
                % Analyze phase for each leg
                for l = 1:nlegs
                    peaks = curr_data.peaksperleg{l};
                    valid_peaks_idx = (peaks >= (refpeak - frame_window)) | ...
                                    (peaks >= refpeak);
                    valid_peaks = peaks(valid_peaks_idx);
                    
                    if isempty(valid_peaks)
                        curr_data.phase_offsets(l, s) = NaN;
                        curr_data.closest_peak(l, s) = NaN;
                        continue;
                    end
                    
                    % Find closest peak to reference
                    [~, min_idx] = min(abs(valid_peaks - refpeak));
                    closest_peak = valid_peaks(min_idx);
                    curr_data.closest_peak(l, s) = closest_peak;
                    
                    % Calculate phase offset
                    phase = mod((closest_peak - refpeak) / period, 1.0);
                    curr_data.phase_offsets(l, s) = phase;
                end
            end
        end
    end
    
    % Update phase_data with computed values
    phase_data(w).closest_peak = curr_data.closest_peak;
    phase_data(w).phase_offsets = curr_data.phase_offsets;
end

end