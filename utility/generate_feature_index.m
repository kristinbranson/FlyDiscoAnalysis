function result = generate_feature_index(feature_names)
% The FlyTracker trk.data array is a 3d double array, fly_count x frame_count x feature_count.
% In order to map from feature names to indices, it's useful to make a scalar struct
% with field names like pos_x, pos_y, etc, each field holding the index into
% the 3d dimension of trk.data.  This function constructs that index struct,
% given as input the trk.names array.

feature_count = numel(feature_names) ;
result = struct() ;
for feature_index = 1 : feature_count ,
  feature_name = feature_names{feature_index} ;
  field_name = strrep(feature_name, ' ', '_') ;
  result.(field_name) = feature_index ;
end

end
