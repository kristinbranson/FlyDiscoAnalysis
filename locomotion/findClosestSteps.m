function [matched_stepdata,matched_idxs] = findClosestSteps(refdata, stepdata, pre_pad_proportion)
% TO DO - does this work on turns? What is duration difference across
% limbs? 
% Find the step onsets closest to the ref limb within the time pretimewindow
% + period of reflimb

% Input
% refdata = set of peaks for reflimb for example that start of swings in a
% walk bout
% stepdata = cell(nlimbs) with the step onsets, peaks, etc 
% reflimb = idx of the reference limb in stepdata
% pre_pad_proportion = porportion of the period before the reference limb
% included in the serach window for closest peaks.

% check to make sure there's more than one step in the walk
% pre_pad_proportion = pi/4;
nlimbs = numel(stepdata);

if ~numel(refdata) < 2

    matched_stepdata = nan(nlimbs,numel(refdata)-1);
    matched_idxs = nan(nlimbs,numel(refdata)-1);
    mperiod = ceil(mean(diff(refdata)));

    for s = 1:numel(refdata)-1
        refstep = refdata(s);
        periods = diff(refdata);
        period = periods(s);
        frame_window = ceil((pre_pad_proportion /(2*pi))*period);

        % % Check if period is reasonable keep??
        % if period >= mperiod * 3
        %     if debug
        %         fprintf('long period frame %d \n',  refpeak);
        %     end
        % else
        % end

        % find closest step onset in other limbs
        for l = 1:nlimbs
            currsteps = stepdata{l};
            valid_step_idx = (currsteps >= (refstep - frame_window)) & ...
                (currsteps < (refstep + period));
            valid_steps = currsteps(valid_step_idx);

            % no close by peaks in the walk bout
            if isempty(valid_steps)
                % matched_stepdata(l, s) = NaN;
                % curr_data.closest_peak(l, s) = NaN;
                continue;
            else

                % Find closest peak to reference
                [~, min_idx] = min(abs(valid_steps - refstep));
                closest_step = valid_steps(min_idx);
                matched_stepdata(l, s) = closest_step;
                matched_idxs(l,s) = find(currsteps == closest_step);
                % check to see if step is used twice
                if ~all(~any(diff(matched_idxs, 1, 2) == 0,2))
                    error('step %d used twice in closest step from leg %d',s,l)
                end
            end

        end


    end
else
    %nan data if only one step in walk ?
end
end



%                 % Calculate phase offset
%                 phase = mod((closest_peak - refpeak) / period, 1.0);
%                 stepdata.phase_offsets(l, s) = phase;
%             end
%         end
%     end
%
% end