function template_file_path = determine_template_or_mask_file_path(specification, rigId, analysis_protocol_folder_path)
% Determines the path to an LED template or mask file, given a specification.
% The specification is something like 'led-image.png', or {'A', 'led-image-A.png',
% 'B', 'led-image-B.png', 'C', 'led-image-C.png', 'D', 'led-image-C.png'}.  In
% the latter case, a different file is used for each rig.  Return value is the
% absolute path to the image file in question.

if ischar(specification) ,
  template_file_path = fullfile(analysis_protocol_folder_path, specification) ;
else
  % In this case, LEDMarkerType should be something like {'A', 'led-image-A.png',
  % 'B', 'led-image-B.png', 'C', 'led-image-C.png', 'D', 'led-image-C.png'}
  rigIdFromRigIndex = specification(1:2:end-1);
  templateFileNameFromRigIndex = specification(2:2:end);
  isMatchFromRigIndex = strcmp(rigId, rigIdFromRigIndex) ;
  matchingRigIndices = find(isMatchFromRigIndex) ;
  matchCount = numel(matchingRigIndices) ;
  if matchCount == 0 ,
    error('specification has no entry for rig %s',rigId) ;
  elseif matchCount == 1,
    matchingRigIndex = matchingRigIndices(1) ;
  else
    error('specification seems to have more than one entry that matches rig %s.  This is no good.', rigId) ;
  end
  template_file_name = templateFileNameFromRigIndex{matchingRigIndex} ;
  template_file_path = fullfile(analysis_protocol_folder_path, template_file_name) ;
end

end
