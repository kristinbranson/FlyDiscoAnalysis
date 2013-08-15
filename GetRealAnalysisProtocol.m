function real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir)

real_analysis_protocol = analysis_protocol;
if isunix,
  try
    while true,
      [status,result] = unix(sprintf('readlink %s',fullfile(settingsdir,real_analysis_protocol)));
      if status,
        break,
      end
      result = strtrim(result);
      if isempty(result),
        break;
      end
      [~,real_analysis_protocol] = myfileparts(result);
    end
  catch ME,
    warning('Error trying to read link location: %s',getReport(ME));
  end
end