function sex = sanitize_metadata_sex(metadata, override_sex, override_gender)
% Determine the proper, modern value for the metadata.sex field, taking
% into account possible overrides, 'sex' (preferred) vs 'gender' (not
% preferred), and various anachronistic ways of encoding the value.
% The returned value is guaranteed to be one of { 'b', 'm', 'f', 'x' }.

% Deal with optional arguments
if ~exist('override_sex', 'var') || isempty(override_sex) ,
  override_sex = '' ;
end
if ~exist('override_gender', 'var') || isempty(override_gender) ,
  override_gender = '' ;
end

% Get the first pass at sex, taking overrides into account.
% We prioritize overrides over metadata, and "sex" over "gender".
if ~isempty(override_sex) ,
  sex0 = override_sex ;
else
  if ~isempty(override_gender),
    sex0 = override_gender;
  else
    if isfield(metadata, 'sex') ,
      sex0 = metadata.sex ;
    elseif isfield(metadata, 'gender') ,
      sex0 = metadata.gender ;
    else
      warningNoTrace('Unable to determine sex from metadata.  Assuming metadata.sex is ''b''') ;
      sex0 = 'b' ;
    end
  end
end

% fix case
sex1 = lower(sex0) ;

% Fix possible metadata error
if strcmp(sex1, 'both') ,
  sex2 = 'b' ;
else
  sex2 = sex1 ;
end
  
% Fix illegal value
valid_sex_values = { 'b', 'm', 'f', 'x' } ;
if any(strcmp(sex2, valid_sex_values)) ,
  sex = sex2 ;
else
  warningNoTrace('metadata.gender is the illegal value ''%s'', setting to ''b''', sex2) ;
  sex = 'b' ;  
end

end  % function
