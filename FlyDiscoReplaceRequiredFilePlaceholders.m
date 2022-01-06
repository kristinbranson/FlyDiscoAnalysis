function required_files_for_stage = ...
        FlyDiscoReplaceRequiredFilePlaceholders(required_files_for_stage_with_placeholders, ...
                                                stage_name, ...
                                                expdir, ...
                                                analysis_protocol_folder_path, ...
                                                dataloc_params, ...
                                                stage_additional_arguments) 
    % Some required files, for some stages, don't have fixed names.  These are
    % typically represented in the required file list for the stage as an all-caps
    % 'placeholder' file name like 'CTRAXRESULTSMOVIE'.  This function replaces
    % these placeholders with the actual file names.  How those file names are
    % determined is custom for each placeholder.

    % If there are no placeholders, we want to return the required file list
    % unaltered.
    required_files_for_stage = required_files_for_stage_with_placeholders ;

    % Handle CTRAXRESULTSMOVIE for stage 'ctraxresultsmovie'
    if strcmp(stage_name, 'ctraxresultsmovie') ,
        i = find(strcmp('CTRAXRESULTSMOVIE',required_files_for_stage),1);
        if ~isempty(i),
            [~,basename] = fileparts(expdir);
            avifilestr = sprintf('%s_%s',dataloc_params.ctraxresultsavifilestr,basename);
            h264file = [avifilestr,'.mp4'];
            required_files_for_stage{i} = h264file;
        end
    end
    
    % Handle PERFRAMEMATFILES for stage 'computeperframefeatures'
    if strcmp(stage_name, 'computeperframefeatures') ,
        if ismember('PERFRAMEMATFILES',required_files_for_stage),
            % read in per-frame fns to compute
            i = find(strcmpi('perframefns',stage_additional_arguments),1);
            if ~isempty(i),
                perframefns = stage_additional_arguments{i+1};
            else
                perframefnsfile = fullfile(analysis_protocol_folder_path, dataloc_params.perframefnsfilestr) ;
                perframefns = importdata(perframefnsfile);
            end
            perframematfiles = cellfun(@(x) fullfile(dataloc_params.perframedir,[x,'.mat']),perframefns,'UniformOutput',false);
            required_files_for_stage = [required_files_for_stage(:)',perframematfiles(:)'];
            required_files_for_stage = setdiff(required_files_for_stage,{'PERFRAMEMATFILES'});
        end
    end
    
    % Handle ANALYSISPLOTS for stage 'plotperframestats'
    if strcmp(stage_name, 'plotperframestats') ,    
        if ismember('ANALYSISPLOTS',required_files_for_stage) ,
            % TODO: Hack alert: we just ignore the 'ANALYSISPLOTS' placeholder for now.  Probably
            % this is not a good long-term solution.
            required_files_for_stage = setdiff(required_files_for_stage, {'ANALYSISPLOTS'}) ;
%             histperframefeaturesfile = fullfile(analysis_protocol_folder_path, dataloc_params.histperframefeaturesfilestr) ;
%             %hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile) ;
%             hist_perframefeatures = ReadStatsPerFrameFeatures2(histperframefeaturesfile);
%             hist_fields = unique({hist_perframefeatures.field}) ;
%             for i = 1:numel(hist_fields) ,
%                 savename = sprintf('hist_%s.png',hist_fields{i});
%                 savename = fullfile(dataloc_params.figdir,savename);
%                 required_files_for_stage{end+1} = savename; %#ok<AGROW>
%             end
%             required_files_for_stage = setdiff(required_files_for_stage, {'ANALYSISPLOTS'}) ;
%             required_files_for_stage{end+1} = fullfile(dataloc_params.figdir, 'stats.png') ;
        end
    end
    
end
