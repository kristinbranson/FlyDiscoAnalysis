function required_files_for_stage = ...
  compute_required_files_for_jaabadetect_stage(stage_name, ...
                                               expdir, ...
                                               analysis_protocol_folder_path, ...
                                               dataloc_params, ...
                                               stage_additional_arguments)  %#ok<INUSD> 

% Determine the output file names for the jaabadetect stage.
% These are typically of the form scores*.mat.

jaabaclassifierparamsfilestrs = fullfile(analysis_protocol_folder_path, dataloc_params.jaabaclassifierparamsfilestrs) ;
raw_jabfiles = read_one_file_name_per_line(jaabaclassifierparamsfilestrs) ;
jabfiles = raw_jabfiles(:)' ;  % want a row vector

function result = score_file_name_from_jab_file_name(jab_file_name)
  Q = loadAnonymous(jab_file_name);
  if isstruct(Q)
    Q = Macguffin(Q);
  end
  Q.modernize(true);  
  raw_score_file_name = Q.file.scorefilename{1} ;
  [~, raw_result] = fileparts(raw_score_file_name) ;
  result = strcat(raw_result, '.mat') ;
end

required_files_for_stage = cellfun(@score_file_name_from_jab_file_name, jabfiles, 'UniformOutput', false) ;

end  % function
