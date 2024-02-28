function result = setfields(s, varargin)
% Utility to set a number of fields of s at once, in a nondestructive way.

if mod(numel(varargin),2) ~= 0 ,
  error('The function %s() only takes an odd number of arguments', mfilename()) ;
end

field_names = varargin(1:2:end) ;
values = varargin(2:2:end) ;

result = s ;
for i = 1 : numel(field_names) ,
  field_name = field_names{i} ;
  value = values{i} ;
  result.(field_name) = value ;
end

