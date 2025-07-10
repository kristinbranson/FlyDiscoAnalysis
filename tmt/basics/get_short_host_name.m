function name = get_short_host_name()
    [ret, name] = system('hostname -s') ;   % Want the short version of the host name, without domain
    if ret == 0 ,
        name = strtrim(name) ;
    else
       if ispc()
          name = getenv('COMPUTERNAME');
       else      
          name = getenv('HOSTNAME');      
       end
       name = strtrim(lower(name));
    end
end
