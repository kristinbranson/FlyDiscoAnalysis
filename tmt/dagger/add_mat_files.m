function add_mat_files(output_file_name, input_file_name_1, input_file_name_2) 
    %error('oops!') ;
    ensure_file_does_not_exist(output_file_name) ;
    input_1 = load_anonymous(input_file_name_1) ;
    input_2 = load_anonymous(input_file_name_2) ;
    output = input_1 + input_2 ;
    save_anonymous(output_file_name, output) ;
end
