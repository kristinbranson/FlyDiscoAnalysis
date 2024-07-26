function result = collect_metadata(experiment_folder_path, metadata_file_name)
  % Collect metadata from various sources, giving preference to more 'modern'
  % sources of metadata.
  
  [metadata_from_experiment_name, did_successfully_extract_metadata_from_experiment_name] = parseExpDir(experiment_folder_path);
  metadata_file_path = fullfile(experiment_folder_path, metadata_file_name) ;
  if exist(metadata_file_path, 'file') ,
    metadata_from_file = ReadMetadataFile(metadata_file_path);
  end  
  if did_successfully_extract_metadata_from_experiment_name ,
    result = structmerge(metadata_from_experiment_name, metadata_from_file) ;  % metadata from file takes precedence
  else
    result = metadata_from_file ;
  end
end



function s = structmerge(varargin)
  % s = structmerge(s1,s2,s3,...)
  % Merge structures into one combined structure.
  % Arguments larther to the right take precedence.

  s = struct();
  for c = 1:numel(varargin)
      stmp = varargin{c};
      assert(isscalar(stmp) && isstruct(stmp),...
          'All input arguments must be scalar structures.');
      tmpflds = fieldnames(stmp);
      for fld = tmpflds(:)'
          fld = fld{1}; %#ok<FXSET>
%           if isfield(s,fld)
%               warning('structmerge:overlappingFields',...
%                   'Overwriting field ''%s''.',fld);
%           end
          s.(fld) = stmp.(fld);
      end
  end
end    
