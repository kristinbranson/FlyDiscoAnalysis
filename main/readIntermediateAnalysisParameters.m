function [intermediate_analysis_parameters, dataloc_params, analysis_protocol_folder_path, datalocparamsfilestr] = ...
  readIntermediateAnalysisParameters(settingsdir, analysis_protocol)
    % Set the default analysis parameters
    default_analysis_parameters_as_list = FlyDiscoPipelineDefaultAnalysisParameters() ;
    default_analysis_parameters = struct_from_name_value_list(default_analysis_parameters_as_list) ;

    % Read in the analysis protocol parameters from the analysis-protocol folder, if
    % it exists
    analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
    analysis_protocol_parameters_file_path = fullfile(analysis_protocol_folder_path, 'analysis-protocol-parameters.txt') ;
    if exist(analysis_protocol_parameters_file_path, 'file') ,
        analysis_protocol_parameters = ReadParams(analysis_protocol_parameters_file_path) ;  % a struct
    else
        analysis_protocol_parameters = struct() ;
    end

    % Read in the dataloc params file, which contains the locations of all the
    % other parameter files
    datalocparamsfilestr = 'dataloc_params.txt' ;
    datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr) ;
    dataloc_params = ReadParams(datalocparamsfile);       
    
    % Backwards-compatibility hack --- If the analysis-protocol folder lacks an
    % analysis-protocol-parameters.txt file, and the indicator_params.txt file
    % exists and says that the experiment is a non-optogenetic experiment, then 
    % turn off the ledonoffdetection stage, since it is not needed, and will in
    % fact error out since FlyDiscoDetectIndicatorLedOnOff() no longer produces an
    % 'empty' output file in this case.
    if is_on_or_force(default_analysis_parameters.doledonoffdetection) ,
      if ~exist(analysis_protocol_parameters_file_path, 'file') ,
        indicatorparamsfile = fullfile(analysis_protocol_folder_path,dataloc_params.indicatorparamsfilestr) ;
        if exist(indicatorparamsfile,'file') ,
          raw_indicator_params = ReadParams(indicatorparamsfile);
          if isfield(raw_indicator_params, 'OptogeneticExp') ,
            if ~logical(raw_indicator_params.OptogeneticExp) ,
              warning(sprintf(['Optogenetic flag is set to false, but the ledonoffdetection stage is enabled by default.\n' ...
                               'You should create an analysis-protocol-parameters.txt in your analysis-protocol file, and turn off this stage.\n' ...
                               'Turning off ledonoffdetection stage.\n'])) ; %#ok<SPWRN> 
              analysis_protocol_parameters.doledonoffdetection = 'off' ;
            end
          end
        end
      end
    end

    % Combine the default parameters with those from the analysis-protocol folder and those in the arguments
    % Precedence is: argument_parameters > analysis-protocol paramters > default parameters
    intermediate_analysis_parameters = merge_structs(default_analysis_parameters, analysis_protocol_parameters) ;  
end