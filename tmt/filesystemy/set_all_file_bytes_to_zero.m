function set_all_file_bytes_to_zero(file_name) 
    % Corrupt a file by changing all bytes to zero.

    [~, ~, byte_count] = simple_dir(file_name) ;
    corrupted_byte_vector = zeros(byte_count, 1, 'uint8') ;
    write_uint8_array_to_file(file_name, corrupted_byte_vector) ;

    %correct_byte_vector = read_file_into_uint8_array(file_name) ;
    %corrupted_byte_vector = correct_byte_vector ;
    %corrupted_byte_vector(:) = uint8(0) ;    
    %byte_count = length(correct_byte_vector)
    %corrupted_byte_vector_check = read_file_into_uint8_array(file_name) ;
    %new_byte_count = length(corrupted_byte_vector_check)     
end
