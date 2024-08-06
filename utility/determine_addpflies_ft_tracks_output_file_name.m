function result = determine_addpflies_ft_tracks_output_file_name(dataloc_params)

if isfield(dataloc_params, 'addpfliesfttracksoutputfilestr') ,
  result = dataloc_params.addpfliesfttracksoutputfilestr ;
else
  result = 'movie-track-with-pflies.mat' ;
end

end
