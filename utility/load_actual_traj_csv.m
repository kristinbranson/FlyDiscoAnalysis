function result = load_actual_traj_csv(filename)

string_from_line_index = read_file_into_cellstring(filename) ;
n = numel(string_from_line_index) ;
if n==0 ,
  result = zeros(0,3) ;
  return
end
result = nan(n-1,3) ;
for i = 2 :n ,
  line = string_from_line_index{i} ;
  tokens = strsplit(line, ',') ;
  if numel(tokens)~=3 ,
    warning('Line %d has %d elements, should have 3.  Skipping.  Row %d in result will be all nans.', i, numel(tokens), i-1) ;
    continue
  end
  result(i-1,1) = str2double(tokens{1}) ;
  result(i-1,2) = str2double(tokens{2}) ;
  token3 = tokens{3} ;
  if strcmp(token3, 'on') ,
    is_on = 1 ;
  elseif strcmp(token3, 'off') ,
    is_on = 0 ;
  else
    warning('In line %d, third token is ''%s'', should be ''on'' or ''off''.  Element (%d,3) in result will be nan.', i, token3, i-1) ;
    is_on = nan ;
  end    
  result(i-1,3) = is_on ;
end
