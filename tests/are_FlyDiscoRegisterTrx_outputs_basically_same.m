function result = are_FlyDiscoRegisterTrx_outputs_basically_same(test_experiment_folder_name, reference_experiment_folder_name)

file_name_from_output_index = { 'registered_trx.mat', 'registrationdata.mat' } ;
output_count = numel(file_name_from_output_index) ;
for output_index = 1 : output_count ,
    output_file_name = file_name_from_output_index{output_index} ;
    test_output_file_path = fullfile(test_experiment_folder_name, output_file_name) ;
    reference_output_file_path = fullfile(reference_experiment_folder_name, output_file_name) ;
    test_output = load('-mat', test_output_file_path) ;
    reference_output = load('-mat', reference_output_file_path) ;
    if strcmp(output_file_name, 'registered_trx.mat') ,
        are_they_same = @are_registered_trx_basically_same ;
    elseif strcmp(output_file_name, 'registrationdata.mat') ,
        are_they_same = @are_registrationdata_basically_same ;
    else
        error('Internal error: Unknown file (%s)', output_file_name) ;    
    end
    this_result = are_they_same(test_output, reference_output) ;
    if ~this_result ,
        result = false ;
        return
    end
end
result = true ;

end
