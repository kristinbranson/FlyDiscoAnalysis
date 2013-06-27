screen_types = {};
for i = 1:numel(tmp5),
  checkfile = fullfile(rootdatadir,tmp5{i}(9:end),'automated_checks_incoming_results.txt');
  metadatafile = fullfile(rootdatadir,tmp5{i}(9:end),'Metadata.xml');
  fprintf('\n%s\n',tmp5{i}(9:end));
  if exist(metadatafile,'file'),
    try
      metadata = ReadMetadataFile(metadatafile);
    catch 
      continue;
    end
    if isfield(metadata,'screen_type'),
      fprintf('screen_type = %s\n',metadata.screen_type);
      screen_types{i} = metadata.screen_type;
    else
      fprintf('no screen type\n');
    end
  else
    fprintf('no metadata\n');
    if ~exist(checkfile,'file'),
      fprintf('no check file\n');
    else
      type(checkfile);
    end
  end
end