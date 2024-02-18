function result = rmfield_quiet(s, field_names_to_delete)
% Annoyingly, rmfield(s, field_names) will error if any of the field_names in
% field_names are not in s.  This function fixes that.

if ~iscell(field_names_to_delete) ,
  field_names_to_delete = { field_names_to_delete } ;
end
original_field_names = fieldnames(s) ;
field_names_to_really_delete = intersect(field_names_to_delete, original_field_names) ;
result = rmfield(s, field_names_to_really_delete) ;
