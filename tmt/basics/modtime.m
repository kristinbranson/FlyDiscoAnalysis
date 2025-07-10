function datenum = modtime(file_name)
    % Returns the modification time of file_name as a Matlab serial date number (days (including fractional days) since some date in the distant past).
    s = dir(file_name) ;
    datenum = s.datenum ;
end
