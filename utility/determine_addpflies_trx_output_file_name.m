function result = determine_addpflies_trx_output_file_name(dataloc_params)

if isfield(dataloc_params, 'addpfliestrxoutputfilestr') ,
  result = dataloc_params.addpfliestrxoutputfilestr ;
else
  result = 'movie_JAABA/trx_with_pflies.mat' ;
end

end
