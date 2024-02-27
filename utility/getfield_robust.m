function result = getfield_robust(s, field_name, fallback)

if isfield(s,field_name) ,
  result = s.(field_name) ;
else
  result = fallback ;
end

