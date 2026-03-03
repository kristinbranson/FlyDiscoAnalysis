function result = flatten_row_cell_array(lst)
% Turn a list (a row cell array) of lists into just a flat list.  Always
% returns a row vector.  Works recursively.  All cell arrays in lst, at any
% depth, should be row cell arrays.
result = cell(1,0) ;
n = numel(lst) ;
for i = 1 : n
  lsti = lst{i} ;
  if iscell(lsti)
    flat_lsti = flatten_row_cell_array(lsti) ;
  else
    flat_lsti = lsti ;
  end
  result = horzcat(result, flat_lsti) ;  %#ok<AGROW>
end
end  % function
