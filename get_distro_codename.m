function result = get_distro_codename()
    if isunix() ,
        [return_code, stdout] = system('lsb_release -i -s') ;
        if return_code==0 ,
            result = strtrim(stdout) ;
        else
            error('Unable to determine operating system---call to lsb_release failed') ;
        end
    else
        if ispc() ,
            result = 'Windows' ;
        elseif ismac() ,
            result = 'macOS' ;
        else
            error('Unable to determine operating system') ;
        end
    end
end
