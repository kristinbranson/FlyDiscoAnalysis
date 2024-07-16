function [out_restrictedSwingStance] = restrictedSwingStance(digital_signal,perfly_limbboutdata,format)
% input options
% format == trajectory is for data in the ctrax trajectory format, 1 = first frame of trajectory
% format == movie is for data in the movie frames format, 1 = first frame of movie
% digtal signal should match the bout data frame of referece

state = {'swing','stance'};


out_restrictedSwingStance = struct;
for fly = 1:numel(perfly_limbboutdata)
    if strcmp(format,'trajectory')
        assert(iscell(digital_signal),'wrong data format')
        curr_digital_signal = digital_signal{fly};
    elseif strcmp(format,'movie')
        assert(isvector(digital_signal), 'wrong data format')
        curr_digital_signal =digital_signal;
    end
    for ileg = 1:numel(perfly_limbboutdata(fly).perlimb)

        for is = 1:(numel(state)) % state = swing or stance
            curr_start_indices = perfly_limbboutdata(fly).perlimb(ileg).(state{is}).start_indices;
            curr_end_indices = perfly_limbboutdata(fly).perlimb(ileg).(state{is}).end_indices;
            [out_start_indices, out_end_indices] = find_bout_overlap(curr_digital_signal, curr_start_indices,curr_end_indices);
            out_restrictedSwingStance(fly).perlimb(ileg).(state{is}).start_indices = out_start_indices;
            out_restrictedSwingStance(fly).perlimb(ileg).(state{is}).end_indices = out_end_indices;

        end

        % stance during digital signal
        %             curr_start_indices = perfly_limbboutdata(fly).perlimb(ileg).start_indices_stance;
        %             curr_end_indices = perfly_limbboutdata(fly).perlimb(ileg).end_indices_stance;
        %
        %             [out_start_indices, out_end_indices] = find_bout_overlap(curr_digital_signal, curr_start_indices,curr_end_indices);
        %             out_restrictedSwingStance(fly).perlimb(ileg).start_indices_stance = out_start_indices;
        %             out_restrictedSwingStance(fly).perlimb(ileg).end_indices_stance = out_end_indices;
        %
        %             % swing during digital signal
        %             curr_start_indices = perfly_limbboutdata(fly).perlimb(ileg).start_indices_swing;
        %             curr_end_indices = perfly_limbboutdata(fly).perlimb(ileg).end_indices_swing;
        %             [out_start_indices, out_end_indices] = find_bout_overlap(curr_digital_signal, curr_start_indices,curr_end_indices);
        %             out_restrictedSwingStance(fly).perlimb(ileg).start_indices_swing = out_start_indices;
        %             out_restrictedSwingStance(fly).perlimb(ileg).end_indices_swing = out_end_indices;
    end
end
