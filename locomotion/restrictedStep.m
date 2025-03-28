function [out_restrictedStep] = restrictedStep(digital_signal,perfly_limbboutdata,format)
% input options
% format == trajectory is for data in the ctrax trajectory format, 1 = first frame of trajectory
% format == movie is for data in the movie frames format, 1 = first frame of movie
% digtal signal should match the bout data frame of referece

% add step data to perfly_limbboutdata - step bouts do NOT line up with
% swing or stance. 


out_restrictedStep = struct;


for fly = 1:numel(perfly_limbboutdata)
    if strcmp(format,'trajectory')
        assert(iscell(digital_signal),'wrong data format')
        curr_digital_signal = digital_signal{fly};
    elseif strcmp(format,'movie')
        assert(isvector(digital_signal), 'wrong data format')
        curr_digital_signal =digital_signal;
    end
    for ileg = 1:numel(perfly_limbboutdata(fly).perlimb)
            % step start = AEP (stance start), step end = next AEP  
            curr_start_indices = perfly_limbboutdata(fly).perlimb(ileg).stance.start_indices(1:end-1);
            curr_end_indices = perfly_limbboutdata(fly).perlimb(ileg).stance.start_indices(2:end);

            [out_start_indices, out_end_indices] = find_bout_overlap(curr_digital_signal, curr_start_indices,curr_end_indices);
            perfly_limbboutdata(fly).perlimb(ileg).step.start_indices = out_start_indices;
            perfly_limbboutdata(fly).perlimb(ileg).step.end_indices = out_end_indices;        
    end
end
out_restrictedStep = perfly_limbboutdata;