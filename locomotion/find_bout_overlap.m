function [out_start_indices, out_end_indices] = find_bout_overlap(digital_signal, in_start_indices,in_end_indices)
% input:
% digital_signal 0,1 array 1 x nframes
% must be paired:
% in_start_indices = start frames of bouts 1 x nbouts
% in_end_indices = end frames of bouts 1 x nbouts

% output: 
% start and end indices for bouts completely within ON of digital signal
  % For example, loop through all stance bouts and return only those where
  % all the frames of the bout are LED on frames.
  out_start_indices = [];
  out_end_indices = [];

  % bouts must be paired
  assert(numel(in_start_indices) == numel(in_end_indices),'bouts indicies must be paired');

  for i = 1:numel(in_start_indices)
      currbout = digital_signal(in_start_indices(i):in_end_indices(i));
      if all(currbout)
          out_start_indices = cat(1,out_start_indices,in_start_indices(i)); 
          out_end_indices = cat(1,out_end_indices,in_end_indices(i));
      end
  end

