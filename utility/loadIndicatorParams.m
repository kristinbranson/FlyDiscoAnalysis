function indicator_params = loadIndicatorParams(file_name) 

raw_indicator_params = ReadParams(file_name);
indicator_params = modernizeIndicatorParams(raw_indicator_params) ;

end
