function real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir)

analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
canonical_analysis_protocol_folder_path = realpath(absolute_filename(analysis_protocol_folder_path)) ;
[~,real_analysis_protocol] = fileparts2(canonical_analysis_protocol_folder_path) ;

end
