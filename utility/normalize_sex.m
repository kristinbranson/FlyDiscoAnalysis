function result = normalize_sex(trx)
% Makes it so that the dimension of trx(i).sex matches the dimensions of
% trx(i).x, trx(i).y, etc.  Errors if it's not clear how to do this.

result = arrayfun(@normalize_sex_for_scalar, trx) ;

end



function result = normalize_sex_for_scalar(trx_scalar)

n = numel(trx_scalar.x) ;
result = trx_scalar ;
result.sex = normalize_sex_for_field(trx_scalar.sex, n) ;
   
end



function result = normalize_sex_for_field(sex, n)

if n==0 ,
  result = cell(1,0) ;
elseif isstringy(sex) ,
  result = repmat({sex}, [1 n]) ;
elseif iscell(sex) ,
  if numel(sex) == n ,
    result = sex ;
  else
    if isscalar(sex) ,      
      result = repmat(sex_string, [1 n]) ;
    else
      error('Unable to normalize sex.') ;
    end
  end
else
  error('Unable to normalize sex.') ;
end  

end  % function
