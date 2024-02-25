function FlyDiscoAutomaticChecksIncoming(expdir, stage, varargin)  %#ok<INUSL>
  % This function's job is to produce the output the files usually
  % named automatic_checks_incoming_results.txt and automatic_checks_incoming_info.mat.
  % If it can't do that, it might error.  If it detects problems with the incoming
  % data that indicate that the pipeline should not proceed, it signals this by 
  % writing 'automated_pf,F' to the text output file, and something similar to the
  % .mat output file, but will still exit without erroring.
  
  version = '0.2.1';
  
  do_continue_pipeline = true;
  msgs = {};
  
  datetime_format = 'yyyymmddTHHMMSS';
  
  [analysis_protocol,settingsdir,datalocparamsfilestr,min_barcode_expdatestr] = ...
    myparse(varargin,...
    'analysis_protocol','current_bubble',...
    'settingsdir', default_settings_folder_path(),...
    'datalocparamsfilestr','dataloc_params.txt',...
    'min_barcode_expdatestr','20110301T000000');
  min_barcode_expdatenum = datenum(min_barcode_expdatestr,datetime_format); %#ok<DATNM> 
  
  %% parameters
  datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
  dataloc_params = ReadParams(datalocparamsfile);
  
  % %%
  % logger = PipelineLogger(expdir,mfilename(),dataloc_params,...
  %   'automaticchecks_incoming_logfilestr',settingsdir,analysis_protocol,...
  %   'logfid',logfid,'debug',DEBUG,'versionstr',version);
  
  try
    paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.automaticchecksincomingparamsfilestr);
    check_params = ReadParams(paramsfile);
    if ~iscell(check_params.control_line_names),
      check_params.control_line_names = {check_params.control_line_names};
    end
    metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
    moviefile = fullfile(expdir,dataloc_params.moviefilestr);
    protocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
    temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr); %#ok<NASGU>
    outfile = fullfile(expdir,dataloc_params.automaticchecksincomingresultsfilestr);
    
    % order matters here: higher up categories have higher priority
    categories = {'flag_aborted_set_to_1',...
      'missing_video',...
      'missing_metadata_file',...
      'missing_metadata_fields',...
      'short_video',...
      'bad_video_timestamps',...
      'missing_capture_files',...
      'flag_redo_set_to_1',...
      'flag_flies_dead_or_damaged',...
      'fliesloaded_time_too_short',...
      'fliesloaded_time_too_long',...
      'shiftflytemp_time_too_long',...
      'no_barcode',...
      'incoming_checks_other'};
    category2idx = struct;
    for i = 1:numel(categories),
      category2idx.(categories{i}) = i;
    end
    iserror = false(1,numel(categories));
    
    %% check for notstarted
    [~,expname] = fileparts(expdir);
    if ~isempty(regexp(expname,'notstarted','once')),
      iserror(category2idx.missing_video) = true;
      do_continue_pipeline = false;
      msgs{end+1} = 'Capture not started';
    end
    
    %% read metadata
    
    ismetadata = exist(metadatafile,'file');
    if ~ismetadata,
      do_continue_pipeline = false;
      msgs{end+1} = 'Missing Metadata.xml file';
      iserror(category2idx.missing_metadata_file) = true;
      metadata = struct;
    else
      try
        metadata = ReadMetadataFile(metadatafile);
      catch %#ok<CTCH>
        msgs{end+1} = 'Error reading Metadata file';
        do_continue_pipeline = false;
      end
    end
    
    
    %% check for metadata fields
    
    required_fns = {'flag_aborted','flag_redo','seconds_fliesloaded','seconds_shiftflytemp',...
      'screen_type','line','cross_barcode'};
    if isfield(check_params,'required_fns'),
      required_fns = check_params.required_fns;
    end
    ismissingfn = ~ismember(required_fns,fieldnames(metadata));
    if any(ismissingfn),
      do_continue_pipeline = false;
      msgs{end+1} = ['Missing required metadata fields:',sprintf(' %s',required_fns{ismissingfn})];
      iserror(category2idx.missing_metadata_fields) = true;
    end
    
    %% check for flags
    
    if isfield(metadata,'flag_aborted') && metadata.flag_aborted ~= 0,
      do_continue_pipeline = false;
      msgs{end+1} = 'Experiment aborted.';
      iserror(category2idx.flag_aborted_set_to_1) = true;
    end
    
    if isfield(metadata,'flag_redo') && metadata.flag_redo ~= 0,
      do_continue_pipeline = false;
      msgs{end+1} = 'Redo flag set to 1.';
      iserror(category2idx.flag_redo_set_to_1) = true;
    end
    
    %% check for dead or damaged flies
    
    if isfield(metadata,'num_flies_damaged') && metadata.num_flies_damaged > 0,
      do_continue_pipeline = false;
      msgs{end+1} = 'Damaged flies > 0.';
      iserror(category2idx.flag_flies_dead_or_damaged) = true;
    end
    
    if isfield(metadata,'num_flies_dead') && metadata.num_flies_dead > 0,
      do_continue_pipeline = false;
      msgs{end+1} = 'Dead flies > 0.';
      iserror(category2idx.flag_flies_dead_or_damaged) = true;
    end
    
    %% check loading time
    
    if isfield(metadata,'seconds_fliesloaded'),
      if metadata.seconds_fliesloaded < check_params.min_seconds_fliesloaded,
        do_continue_pipeline = false;
        msgs{end+1} = sprintf('Load time = %f < %f seconds.',metadata.seconds_fliesloaded,check_params.min_seconds_fliesloaded);
        iserror(category2idx.fliesloaded_time_too_short) = true;
      end
      if metadata.seconds_fliesloaded > check_params.max_seconds_fliesloaded,
        do_continue_pipeline = false;
        msgs{end+1} = sprintf('Load time = %f > %f seconds.',metadata.seconds_fliesloaded,check_params.max_seconds_fliesloaded);
        iserror(category2idx.fliesloaded_time_too_long) = true;
      end
    end
    
    %% check shiftflytemp time
    % if isfield(metadata,'seconds_shiftlytemp') && ...
    %     metadata.seconds_shiftflytemp > check_params.max_seconds_shiftflytemp,
    %   success = false;
    %   msgs{end+1} = sprintf('Shift fly temp time = %f > %f seconds.',metadata.seconds_shiftflytemp,check_params.max_seconds_shiftflytemp);
    %   iserror(category2idx.shiftflytemp_time_too_long) = true;
    % end
    
    %% check for primary screen: all other checks are only for screen data
    
    if ~isfield(metadata,'screen_type'),
      % should already be captured in missing metadata fields
      %   success = false;
      %   msgs{end+1} = 'screen_type not stored in Metadata file';
      isscreen = false;
    else
      isscreen = strcmpi(metadata.screen_type,'non_olympiad_branson_cx');
    end
    
    %% check barcode
    
    if isfield(metadata,'exp_datetime') && ~isempty(metadata.exp_datetime) ,
      exp_datenum = datenum(metadata.exp_datetime,datetime_format); %#ok<DATNM> 
    else
      exp_datenum = nan;
    end
    
    
    %% check video length
    
    headerinfo = [];
    
    registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
    if ~exist(registrationparamsfile,'file'),
      error('Registration params file %s does not exist',registrationparamsfile);
    end
    registration_params = ReadParams(registrationparamsfile);
    
    indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
    if exist(indicatorparamsfile,'file') ,
      indicator_params = ReadParams(indicatorparamsfile) ;
      if ~isfield(indicator_params, 'OptogeneticExp') ,
        error('No optogenetics flag in indicator params')
      end
    else
      error('Indicator params file %s does not exist',indicatorparamsfile);
    end

    if exist(moviefile,'file')
      try
        headerinfo = ufmf_read_header(moviefile);
        % if optogenetic experiment
        if indicator_params.OptogeneticExp
          % if uses mediandt = unreliable timestamps
          if registration_params.usemediandt
            if ~isfield(headerinfo,'nframes'),
              error('No field nframes in UFMF header');
            end
            nframes = headerinfo.nframes;
            movielength = nframes/check_params.frame_rate;
          else % use timestamps
            if ~isfield(headerinfo,'timestamps'),
              error('No field timestamps in UFMF header');
            end
            movielength = headerinfo.timestamps(end);
          end
          if exist(protocolfile,'file')  % protocol file is checked for below not writing error here if doesn't exist
            protocolinfo = load(protocolfile);
            if ~isfield(protocolinfo.protocol,'duration')
              error('No duration field in protocol file')
            end
            movielength_expected = sum(protocolinfo.protocol.duration)/1000;
            if isfield(check_params,'movie_length_delta_seconds')
              movielength_expected = movielength_expected  - check_params.movie_length_delta_seconds;
            end
            if movielength < movielength_expected
              do_continue_pipeline = false;
              msgs{end+1} = sprintf('Video duration is %d sec < %d sec the protocol duration.',movielength,movielength_expected);
              iserror(category2idx.short_video) = true;
            end
          end
        elseif isfield(check_params,'min_movie_length_seconds')
          if ~isfield(headerinfo,'timestamps'),
            error('No field timestamps in UFMF header');
          end
          movielength = headerinfo.timestamps(end);
          
          movielength_expected = check_params.min_movie_length_seconds;
          
          if movielength < movielength_expected
            do_continue_pipeline = false;
            msgs{end+1} = sprintf('Video duration is %d sec < %d sec min_movie_length.',movielength,movielength_expected);
            iserror(category2idx.short_video) = true;
          end
        elseif isfield(check_params,'min_ufmf_diagnostics_summary_nframes')
          if ~isfield(headerinfo,'nframes'),
            error('No field nframes in UFMF header');
          end
          nframes = headerinfo.nframes;
          
          min_nframes = check_params.min_ufmf_diagnostics_summary_nframes;
          
          if nframes < min_nframes
            do_continue_pipeline = false;
            msgs{end+1} = sprintf('Video nframes is %d < %d the minium frame length.',nframes, min_nframes);
            iserror(category2idx.short_video) = true;
          end
        end
      catch ME,
        do_continue_pipeline = false;
        msgs{end+1} = sprintf('Error reading data from movie header: %s',getReport(ME,'extended','hyperlinks','off'));
        iserror(category2idx.incoming_checks_other) = true;
      end
    end
    
    %% check for bad timestamps
    
    if ~isempty(headerinfo),
      
      if ~isfield(headerinfo,'timestamps'),
        msgs{end+1} = 'No field timestamps in UFMF header';
        do_continue_pipeline = false;
        iserror(category2idx.incoming_checks_other) = true;
      else
        dt = diff(headerinfo.timestamps);
        nneg = nnz(dt <= 0);
        if nneg > 0,
          do_continue_pipeline = false;
          msgs{end+1} = sprintf('%d frames in movie with non-increasing timestamps',nneg);
          iserror(category2idx.bad_video_timestamps) = true;
        end
      end
      
    end
    
    %%
    
    if isscreen,
      
      %exp_datenum = datenum(metadata.exp_datetime,datetime_format);
      
      if (~ismember(metadata.line,check_params.control_line_names) || ...
          (exp_datenum >= min_barcode_expdatenum)) && ...
          metadata.cross_barcode < 0,
        do_continue_pipeline = false;
        msgs{end+1} = 'Barcode = -1 and line_name indicates not control';
        iserror(category2idx.no_barcode) = true;
      end
      
      %% check temperature
      %
      % if ~exist(temperaturefile,'file'),
      %   warning('Temperature file %s does not exist',temperaturefile);
      % else
      %   try
      %     tempdata = importdata(temperaturefile,',');
      %   catch ME,
      %     warning('Error importing temperature stream: %s',getReport(ME));
      %   end
      %   if isempty(tempdata),
      %     warning('No temperature readings recorded.');
      %     temp = [];
      %   elseif size(tempdata,2) < 2,
      %     warning('Temperature data could not be read');
      %     temp = [];
      %   else
      %     temp = tempdata(:,2);
      %     if isempty(temp),
      %       warning('No temperature readings recorded.');
      %     else
      %       if max(temp) > check_params.max_temp,
      %         success = false;
      %         msgs{end+1} = sprintf('Max temperature = %f > %f.',max(temp),check_params.max_temp);
      %       end
      %       if numel(temp) < 2,
      %         warning('Only one temperature recorded.');
      %       else
      %         if max(temp) - min(temp) > check_params.max_tempdiff,
      %           success = false;
      %           msgs{end+1} = sprintf('Temperature change = %f > %f.',max(temp) - min(temp),check_params.max_tempdiff);
      %         end
      %       end
      %     end
      %   end
      % end
      
    end
    
    %% check for missing files
    
    fn = dataloc_params.moviefilestr;
    isfile = exist(fullfile(expdir,fn),'file');
    if ~isfile,
      iserror(category2idx.missing_video) = true;
    end
    if isfield(check_params,'required_files_mindatestr'),
      mindatenum = nan(size(check_params.required_files_mindatestr));
      for i = 1:numel(check_params.required_files_mindatestr),
        if ~isempty(check_params.required_files_mindatestr{i}),
          mindatenum(i) = datenum(check_params.required_files_mindatestr{i},datetime_format); %#ok<DATNM> 
        end
      end
    else
      mindatenum = [];
    end
    
    for i = 1:numel(check_params.required_files),
      if i <= numel(mindatenum) && ~isnan(mindatenum(i)) && ...
          (isnan(exp_datenum) || exp_datenum < mindatenum(i)),
        continue;
      end
      fn = check_params.required_files{i};
      if any(fn == '*'),
        isfile = ~isempty(dir(fullfile(expdir,fn)));
      else
        isfile = exist(fullfile(expdir,fn),'file');
      end
      if ~isfile,
        msgs{end+1} = sprintf('Missing file %s',fn); %#ok<AGROW>
        do_continue_pipeline = false;
        iserror(category2idx.missing_capture_files) = true;
      end
    end
    
    %% check for missing desired files
    
    if isfield(check_params,'desired_files'),
      for i = 1:numel(check_params.desired_files),
        fn = check_params.desired_files{i};
        if any(fn == '*'),
          isfile = ~isempty(dir(fullfile(expdir,fn)));
        else
          isfile = exist(fullfile(expdir,fn),'file');
        end
        if ~isfile,
          msgs{end+1} = sprintf('Missing desired file %s',fn); %#ok<AGROW>
        end
      end
    end
    
    %% output results to text file
    % First just rite text-file contents to a string
    if do_continue_pipeline ,
      text_output_file_as_string = sprintf('automated_pf,U\n');
    else
      text_output_file_as_string = sprintf('automated_pf,F\n');
      i = find(iserror,1);
      if isempty(i),
        s = 'incoming_checks_other';
      else
        s = categories{i};
      end
      text_output_file_as_string = horzcat(text_output_file_as_string,sprintf('automated_pf_category,%s\n',s)) ;
    end
    if ~isempty(msgs),
      text_output_file_as_string = horzcat(text_output_file_as_string,'notes_curation,') ;
      s = sprintf('%s\\n',msgs{:});
      text_output_file_as_string = horzcat(text_output_file_as_string,s) ;
    end
    % Now write the text-file-as-string to the file
    if exist(outfile,'file'),
      try
        delete(outfile);
      catch ME,
        warning('FlyBubbleAutomaticChecksIncoming:output',...
          'Could not delete pre-existing file %s:\n %s',outfile,getReport(ME));
      end
    end
    try
      write_string_to_text_file(outfile, text_output_file_as_string) ;
    catch ME,
      warning('FlyBubbleAutomaticChecksIncoming:output',...
              'Could not write file %s:\n %s',outfile,getReport(ME));
    end        
  catch ME,
    do_continue_pipeline = false;
    msgs = {getReport(ME)};
  end
  
  %% save info to mat file
  
  filename = fullfile(expdir,dataloc_params.automaticchecksincominginfomatfilestr);
  fprintf('Saving debug info to file %s...\n',filename);
  
  %aciinfo = runInfo;
  aciinfo = struct() ;
  aciinfo.paramsfile = paramsfile;
  aciinfo.check_params = check_params;
  aciinfo.version = version;
  aciinfo.iserror = iserror;
  aciinfo.categories = categories;
  aciinfo.msgs = msgs;
  aciinfo.success = do_continue_pipeline;
  
  if exist(filename,'file'),
    try %#ok<TRYNC>
      delete(filename);
    end
  end
  try
    save(filename,'-struct','aciinfo');
  catch ME,
    warning('FlyBubbleAutomaticChecksIncoming:save',...
      'Could not save information to file %s: %s',filename,getReport(ME));
  end  
  
  % Print success/failure and messages to stdout
  fprintf('success = %d\n',do_continue_pipeline);
  if isempty(msgs),
    fprintf('No error or warning messages.\n');
  else
    fprintf('Warning/error messages:\n');
    fprintf('%s\n',msgs{:});
  end
  %logger.close();
  
%   % If not successful, throw an error to stop the pipeline
%   if ~success ,
%     flydisco_pipeline_error(stage, msgs) ;
%   end
end

