function result = is_on_or_force(on_off_force)
  result = isequal(on_off_force, 'on') || isequal(on_off_force, 'force') ;
end
