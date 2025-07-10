function corrupt(file_name) 
    % Corrupt a file by changing a random byte to a random value different from
    % the original value
    correct_byte_vector = read_file_into_uint8_array(file_name) ;
    len = length(correct_byte_vector) ;
    i = randi(len, 1) ;
    correct_value = correct_byte_vector(i) ;
    corrupt_value = uint8(randi(256,1)-1) ;
    while corrupt_value == correct_value ,
        corrupt_value = uint8(randi(256,1)-1) ;
    end
    corrupted_byte_vector = correct_byte_vector ;
    corrupted_byte_vector(i) = corrupt_value ;
    write_uint8_array_to_file(file_name, corrupted_byte_vector) ;
end
