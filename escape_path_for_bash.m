function result = escape_path_for_bash(path)
    % This will work on any string, really, not just a string containing a path
    sq = '''' ;  % a single quote
    sq_bs_sq_sq = '''\''''' ; % single quote, backslash, single quote, single quote
    path_escaped = strrep(path, sq, sq_bs_sq_sq) ;  % replace ' with '\''
    result = horzcat(sq, path_escaped, sq) ;  % Surround with single quotes to handle all special chars besides single quote
end
