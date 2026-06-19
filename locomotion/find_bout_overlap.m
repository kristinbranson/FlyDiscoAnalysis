function [out_start_indices, out_end_indices] = find_bout_overlap(digital_signal, in_start_indices,in_end_indices,end_exclusive)
% input:
% digital_signal 0,1 array 1 x nframes
% must be paired:
% in_start_indices = start frames of bouts 1 x nbouts
% in_end_indices = end frames of bouts 1 x nbouts
% end_exclusive (optional, default false):
%   true  - in_end_indices is the EXCLUSIVE end (first frame AFTER the bout,
%           as for swing/stance from detect_bouts); overlap tested over
%           start:end-1 (the bout's real frames).
%   false - in_end_indices is the last real frame of the bout (as for step,
%           where end = next touchdown); overlap tested over start:end.
%   The RETURNED end indices are always the original in_end_indices, unchanged.

% output:
% start and end indices for bouts completely within ON of digital signal
  % For example, loop through all stance bouts and return only those where
  % all the frames of the bout are LED on frames.
  if nargin < 4 || isempty(end_exclusive)
      end_exclusive = false;
  end
  out_start_indices = [];
  out_end_indices = [];

  % bouts must be paired
  assert(numel(in_start_indices) == numel(in_end_indices),'bouts indicies must be paired');

  for i = 1:numel(in_start_indices)
      last_frame = in_end_indices(i);
      if end_exclusive
          last_frame = last_frame - 1;   % swing/stance: drop the first-frame-after
      end
      currbout = digital_signal(in_start_indices(i):last_frame);
      if all(currbout)
          out_start_indices = cat(1,out_start_indices,in_start_indices(i));
          out_end_indices = cat(1,out_end_indices,in_end_indices(i));   % original, unchanged
      end
  end

