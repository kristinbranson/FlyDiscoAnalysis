function result = FlyDiscoComputeRequiredFileNamesForAllStages(analysis_protocol_folder_path, dataloc_params, expdir, do_run)

% Determine what all the file names are that should be produced by each
% pipeline stage.  Returns a struct s such that s.(stage_name) is a cell
% string of the required file names for each stage.  The field names of s are
% those returned by FlyDiscoStageNames().  All returned file names are relative to the
% experiment folder.

% This is sometimes needed to determine what files to expect
metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
metadata = ReadMetadataFile(metadatafile);

name_from_stage_index = FlyDiscoStageNames() ;
stage_count = numel(name_from_stage_index) ;

result = struct() ;
for stage_index = 1 : stage_count ,
  stage_name = name_from_stage_index{stage_index} ;
  if strcmp(do_run.(stage_name), 'off') ,
    % Skip stages that are turned off, b/c the relevant dataloc_params fields
    % might be missing.
    continue
  end
  if strcmp(stage_name, 'automaticchecksincoming') ,
    required_file_names = { dataloc_params.automaticchecksincomingresultsfilestr, ...
                            dataloc_params.automaticchecksincominginfomatfilestr } ;
  elseif strcmp(stage_name, 'flytracking') ,
    required_file_names = { dataloc_params.flytrackercalibrationstr, ...
                            dataloc_params.flytrackertrackstr, ...
                            dataloc_params.flytrackerbgstr, ...
                            dataloc_params.ctraxfilestr, ...
                            'movie-feat.mat', ...
                            'movie-params.mat'} ;
  elseif strcmp(stage_name, 'addpflies') ,
    addpflies_trx_output_file_name = determine_addpflies_trx_output_file_name(dataloc_params) ;
    addpflies_ft_tracks_output_file_name = determine_addpflies_ft_tracks_output_file_name(dataloc_params) ;
    required_file_names = { addpflies_trx_output_file_name, ...
                            addpflies_ft_tracks_output_file_name } ;
  elseif strcmp(stage_name, 'registration') ,
    required_file_names = { dataloc_params.trxfilestr, ...
                            dataloc_params.registrationmatfilestr, ...
                            dataloc_params.registrationtxtfilestr, ...
                            dataloc_params.registrationimagefilestr } ;    
  elseif strcmp(stage_name, 'ledonoffdetection') ,
    required_file_names = { dataloc_params.indicatordatafilestr } ;    
  elseif strcmp(stage_name, 'sexclassification') ,
    sex = sanitize_metadata_sex(metadata) ;
    if strcmp(sex,'b') ,
      required_file_names = { dataloc_params.sexclassifierdiagnosticsfilestr dataloc_params.sexclassifiermatfilestr } ;
    else
      required_file_names = { dataloc_params.sexclassifierdiagnosticsfilestr } ;
    end        
    % N.B.: This stage also modifies registered_trx.mat
  elseif strcmp(stage_name, 'computeperframefeatures') ,
    perframe_dir_name = dataloc_params.perframedir ;
    perframefnsfilename = fullfile(analysis_protocol_folder_path, dataloc_params.perframefnsfilestr) ;
    perframefns = importdata(perframefnsfilename);    
    perframematfiles = cellfun(@(x) fullfile(perframe_dir_name,[x,'.mat']),perframefns,'UniformOutput',false);
    required_file_names = horzcat({perframe_dir_name}, perframematfiles(:)') ;
  elseif strcmp(stage_name, 'computehoghofperframefeatures') ,
    required_file_names = { 'hs_sup_01_01_1.mat', 'hf_01_01_1.mat' } ;
  elseif strcmp(stage_name, 'jaabadetect') ,
    required_file_names = compute_required_files_for_jaabadetect_stage(stage_name, ...
                                                                       expdir, ...
                                                                       analysis_protocol_folder_path, ...
                                                                       dataloc_params) ;
  elseif strcmp(stage_name, 'locomotionmetrics') ,
    required_file_names = { 'perframe/nfeet_ground.mat', ...
                            dataloc_params.locomotionmetricsswingstanceboutstatsfilestr } ;
  elseif strcmp(stage_name, 'computeperframestats') ,
    required_file_names = { dataloc_params.statsperframetxtfilestr, ...
                            dataloc_params.statsperframematfilestr, ...
                            dataloc_params.histperframetxtfilestr, ...
                            dataloc_params.histperframematfilestr} ;
  elseif strcmp(stage_name, 'plotperframestats') ,
    required_file_names = { dataloc_params.figdir, dataloc_params.statswebpagefilestr } ;    
    histperframefeaturesfile = fullfile(analysis_protocol_folder_path, dataloc_params.histperframefeaturesfilestr) ;
    %hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile) ;
    hist_perframefeatures = ReadStatsPerFrameFeatures2(histperframefeaturesfile);
    hist_fields = unique({hist_perframefeatures.field}) ;
    for i = 1:numel(hist_fields) ,
      savename = sprintf('hist_%s.png',hist_fields{i});
      savename = fullfile(dataloc_params.figdir,savename);
      required_file_names{end+1} = savename; %#ok<AGROW>
    end
  elseif strcmp(stage_name, 'makectraxresultsmovie') ,
    [~,basename] = fileparts(expdir) ;
    avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr, basename) ;
    h264filename = horzcat(avifilestr, '.mp4') ;
    required_file_names = { h264filename } ;
  elseif strcmp(stage_name, 'apt') ,
    required_file_names = { dataloc_params.apttrkfilestr, ...
                            dataloc_params.aptlogfilestr, ...
                            dataloc_params.apterrfilestr, ...
                            dataloc_params.aptbsublogfilestr ...
                            } ;
  elseif strcmp(stage_name, 'makeaptresultsmovie') ,
    [~,basename] = fileparts(expdir) ;
    mp4filestr = sprintf('%s_%s', dataloc_params.aptresultsavefilestr, basename) ;
    apth264file = horzcat(mp4filestr, '.mp4') ;
    required_file_names = { apth264file } ;
  elseif strcmp(stage_name, 'automaticcheckscomplete') ,
    required_file_names = { dataloc_params.automaticcheckscompleteresultsfilestr, ...
                            dataloc_params.automaticcheckscompleteinfomatfilestr } ;
  elseif strcmp(stage_name, 'cleanup') ,
    required_file_names = cell(1,0) ;
  else
    error('Unknown stage name %s', stage_name) ;
  end
  
  % Store the required files in the struct
  result.(stage_name) = required_file_names ;
end

end
