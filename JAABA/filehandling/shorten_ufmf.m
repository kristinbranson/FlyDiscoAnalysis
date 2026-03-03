function shorten_ufmf(input_file_name, output_file_name, desired_frame_count, first_input_frame)
  % SHORTEN_UFMF  Extract a range of frames from a .ufmf file.
  %   shorten_ufmf(input, output, count) extracts the first count frames.
  %   shorten_ufmf(input, output, count, first_input_frame) extracts count
  %   frames starting from frame first_input_frame (1-indexed).

  if nargin < 4 ,
    first_input_frame = 1 ;
  end

  % read in the header
  input_header = ufmf_read_header(input_file_name) ;
  input_frame_count = input_header.nframes ;

  % find the end of the header
  headerstart = 0;
  % loc of something after the header
  headerend = min(min(input_header.mean2file),min(input_header.frame2file));

  % read in the header
  fseek(input_header.fid,headerstart,'bof');
  input_header_as_bytes = fread(input_header.fid,headerend-headerstart);

  available_frame_count = input_frame_count - first_input_frame + 1 ;
  if available_frame_count <= 0 ,
    error('first_input_frame (%d) is > the number of frames in the input (%d)', first_input_frame, input_frame_count) ;
  end
  output_frame_count = min(available_frame_count, desired_frame_count) ;

  last_input_frame = first_input_frame + output_frame_count - 1 ;

  if exist(output_file_name,'file') ,
    delete(output_file_name) ;
  end

  output_fid = fopen(output_file_name, 'wb' , 'ieee-le') ;
  if output_fid < 0 ,
    error('Unable to open output file %s', output_file_name) ;
  end
  cleaner = onCleanup(@()(fclose(output_fid))) ;

  % copy the header to the output
  fwrite(output_fid, input_header_as_bytes) ;

  lastmean = nan;

  index = struct;
  index.frame = struct;
  index.frame.loc = int64([]);
  index.frame.timestamp = input_header.timestamps(first_input_frame:last_input_frame);
  index.keyframe = struct;
  index.keyframe.mean = struct() ;  % this level seems superfluous, but ufmf_read_header() seems to want it
  index.keyframe.mean.loc = int64([]);
  index.keyframe.mean.timestamp = [];


  % write the frames
  for input_frame_index = first_input_frame:last_input_frame ,
    % write a mean if this frame uses a different mean than the last
    if input_header.frame2mean(input_frame_index) ~= lastmean,
      meanstart = input_header.mean2file(input_header.frame2mean(input_frame_index));
      meanend = min([min(input_header.mean2file(input_header.mean2file>meanstart)),...
                     min(input_header.frame2file(input_header.frame2file>meanstart)),...
                     input_header.indexloc]);
      assert(~isempty(meanend) && ~isnan(meanend));
      index.keyframe.mean.loc(end+1) = ftell(output_fid);
      index.keyframe.mean.timestamp(end+1) = input_header.mean_timestamps(input_header.frame2mean(input_frame_index));
      fseek(input_header.fid,meanstart,'bof');
      fwrite(output_fid,fread(input_header.fid,meanend-meanstart));
      lastmean = input_header.frame2mean(input_frame_index);
    end

    framestart = input_header.frame2file(input_frame_index);
    frameend = min([min(input_header.mean2file(input_header.mean2file>framestart)),...
                    min(input_header.frame2file(input_header.frame2file>framestart)),...
                    input_header.indexloc]);
    assert(~isempty(frameend) && ~isnan(frameend));
    fseek(input_header.fid,framestart,'bof');
    index.frame.loc(end+1) = ftell(output_fid);
    fwrite(output_fid,fread(input_header.fid,frameend-framestart));
  end

  % write the index
  index_offset = ftell(output_fid) ;
  ufmf_write_struct(output_fid, index) ;

  % overwrite the index location in the header
  offset = input_header.indexlocloc ;  % The offset where we should write the index offset
  fseek(output_fid, offset, 'bof') ;
  fwrite(output_fid, index_offset, 'uint64') ;
end
