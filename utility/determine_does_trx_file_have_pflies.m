function result = determine_does_trx_file_have_pflies(file_name) 
% Test whether a trx file has any pflies in it.

trx_wrapper = load('-mat', file_name) ;
trx = trx_wrapper.trx ;
if isfield(trx, 'is_pfly') ,
  is_pfly = [trx.is_pfly] ;
  result =  any(is_pfly) ;
else
  result = false ;
end
