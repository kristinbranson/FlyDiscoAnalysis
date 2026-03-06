function [trx, registration_data] = ...
  performTemporalTruncation(dotemporaltruncation, trx, registration_data, newid2oldid, fns_notperframe, recordLengthIdeal, moviefile)

% Truncate end of movie based on value from registration params.

if dotemporaltruncation ,
  ndelete = 0;
  
  headerInfo_local = ufmf_read_header(moviefile);
  if ~isfield(headerInfo_local,'timestamps'),
    error('No field timestamps in UFMF header');
  end
  timestamps_header = headerInfo_local.timestamps;    %??? couldn't figure out where timestamps come from in load_tracks for movie_JAABA/trx.mat
    
  % how long is the video
  recordLengthCurr = timestamps_header(end);
    
  % how much time should we crop from the end?
  if recordLengthCurr < recordLengthIdeal,
    warning('Cropped video is %f seconds long, shorter than ideal length %f seconds.',recordLengthCurr,recordLengthIdeal);
    i1 = numel(timestamps_header);
  else
    i1 = find(timestamps_header >= recordLengthIdeal,1);
    if isempty(i1),
      warning('No timestamps occur after recordLengthIdeal = %f. Cannot crop end.',recordLengthIdeal); %??? why warn and not error? 
      i1 = numel(timestamps_header);
    else
      %??? couldn't figure out what this is checking for
      if i1 > 1 && ...
          (recordLengthIdeal - timestamps_header(i1-1)) < ...
          (timestamps_header(i1) - recordLengthIdeal), 
        i1 = i1 - 1;
      end
    end
  end

  %not cropping start
  i0 = 1;
  registration_data.seconds_crop_start = 0;
  registration_data.start_frame = i0;
  registration_data.seconds_crop_end = timestamps_header(end)-timestamps_header(i1);
  registration_data.end_frame = i1;
  
  fns = setdiff(fieldnames(trx),fns_notperframe);
  isperframe = true(1,numel(fns));
  nperfn = nan(1,numel(fns));
  if ~isempty(trx),    
    for j = 1:numel(fns),
      fn = fns{j};
      for i = 1:numel(trx),
        if ~isnumeric(trx(i).(fn)),
          isperframe(j) = false;
          break;
        end
        % check if fn has perframe data 
        if i == 1,
          ncurr = trx(i).nframes - numel(trx(i).(fn)); 
        else
          if trx(i).nframes - numel(trx(i).(fn)) ~= ncurr,
            isperframe(j) = false;
            break;
          end
        end
      end
      if all([trx.nframes]) == trx(1).nframes && ...
          trx(1).nframes > 1 && numel(trx(1).(fn)) == 1,
        isperframe(j) = false;
      end
      if isperframe(j),
        nperfn(j) = ncurr;
      end
    end
    
    nperfn = nperfn(isperframe);
    fns = fns(isperframe);
    ncropright = ceil(nperfn/2);
    ncropleft = nperfn - ncropright;
    
    trxdelete = false(1,numel(trx));
    for i = 1:numel(trx),
      if trx(i).firstframe > i1,
        trxdelete(i) = true;
        continue;
      end
      
      trx(i).nframes = min(i1,trx(i).endframe)-max(i0,trx(i).firstframe)+1;
      
      if trx(i).firstframe < i0,
        off = i0 - trx(i).firstframe;
        for j = 1:numel(fns),
          fn = fns{j};
          trx(i).(fn) = trx(i).(fn)(off+1+ncropleft(j):end);
        end
        trx(i).firstframe = i0;
      end
      
      if trx(i).endframe > i1,
        for j = 1:numel(fns),
          fn = fns{j};
          trx(i).(fn) = trx(i).(fn)(1:trx(i).nframes-nperfn(j));
        end
        trx(i).endframe = i1;
      end
      
      trx(i).off = -trx(i).firstframe + 1;      
    end
    trx(trxdelete) = []; 
    newid2oldid(trxdelete) = [];
    ndelete = nnz(trxdelete);    
  end
end

% Need to stuff these things into registration_data whether or not we do
% temporal truncation
registration_data.newid2oldid = newid2oldid;
registration_data.flytracker_nidsnew = numel(trx);

% Print an informative message
if dotemporaltruncation ,
  fprintf('Applied temporal truncation. Data cropped to %f seconds. Deleted %d trajectories.\n',timestamps_header(i1),ndelete);
else
  fprintf('NOT applying temporal truncation.\n')
end

end  % function
