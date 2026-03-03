function result = is_m_file(file_name) 
    if length(file_name)<2 ,
        result = false ;
    else
        result = strcmp(file_name(end-1:end), '.m') ;
    end
end
