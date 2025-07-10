function increment_mat_file(output_file_name, input_file_name) 
    %error('oops!') ;
    ensure_file_does_not_exist(output_file_name) ;
    input = load_anonymous(input_file_name) ;
    output = input + 1 ;
    save_anonymous(output_file_name, output) ;
end
