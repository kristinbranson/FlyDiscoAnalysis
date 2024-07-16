function [durations_frames,durations_time] = computeBoutDurations(start_indices, end_indicies,timestamps)
% output framse and time in milliseconds
durations_frames = (end_indicies-start_indices)';
durations_time = timestamps(end_indicies)-timestamps(start_indices);
durations_time = durations_time.*1000;