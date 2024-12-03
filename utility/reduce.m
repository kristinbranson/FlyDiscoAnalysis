function result = reduce(fn, lst, initial_value)
% A reduce function, like in Lisp.  Applies the two-argument function fn to
% initial_value and lst(1), yielding result_1.  Then applies the two-argument function fn to
% result_1 and lst(2), yielding result_2.  Continues until all elements of lst
% have been 'consumed'.  Final result is result_n, n the length of lst.  For
% an empty list, returns initial_value.  If lst is a cell array, uses lst{i}
% to access the i'th element of lst.  Otherwise, uses lst(i).

n = numel(lst) ;
result = initial_value ;
if iscell(lst) ,
  for i = 1:n ,
    result = feval(fn, result, lst{i}) ;
  end
else
  for i = 1:n ,
    result = feval(fn, result, lst(i)) ;
  end
end  

end
