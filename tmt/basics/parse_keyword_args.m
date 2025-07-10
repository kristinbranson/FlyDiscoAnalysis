function varargout = parse_keyword_args(key_value_pairs, varargin)
    % Parse an arbitrary number of keyword arguments.
    % key_value_pairs should be a list of key-value pairs, as supplied when
    % calling a function.  varargin should contain a list of key-default pairs.
    % On return, varargout will contain the values from key_value_pairs, but in
    % the order indicated by varargin.  Any keys missing from key_value_pairs will
    % get the default specified in the even-indexed elements of varargin.
    %
    % This is basically a copy of Kristin Branson's myparse(), which includes
    % optimization(s) by Ben Arthur.

    key_default_pairs = varargin ;  % wasteful, but who cares?
    key_default_pairs_length = length(key_default_pairs) ;
    if mod(key_default_pairs_length, 2) ~= 0 ,
        error('The list of key-default pairs must have an even number of elements') ;
    end
    valid_key_count = round(key_default_pairs_length/2) ;
    key_value_pairs_length = length(key_value_pairs) ;
    if mod(key_value_pairs_length, 2) ~= 0 ,
        error('The list of key-value pairs must have an even number of elements') ;
    end
    key_count = round(key_value_pairs_length/2) ;    
    if nargout ~= valid_key_count ,
        error('The number of output arguments must be equal to the number of valid keys') ;
    end

    % Break out the valid keys, and the defauls
    key_from_valid_key_index = key_default_pairs(1:2:end) ;
    default_from_valid_key_index = key_default_pairs(2:2:end) ;

    % Set all outputs to the defaults
    varargout = default_from_valid_key_index ;

    % Break out the keys, values
    key_from_key_index = key_value_pairs(1:2:end) ;
    value_from_key_index = key_value_pairs(2:2:end) ;

    % Lookup each key in the valid-key list
    for key_index = 1 : key_count ,
        % Get the key, and the value
        key = key_from_key_index{key_index} ;
        value = value_from_key_index{key_index} ;

        % Check that the key is a string
        if ~ischar(key),
            error('Input %d is not a string', 2*key_index+1);
        end
        
        % Find the matching elements in the valid key list
        matching_valid_key_indices = find(strcmp(key, key_from_valid_key_index)) ;

        % Deal with possible misses, multiple matches
        match_count = length(matching_valid_key_indices) ;
        if match_count == 0 ,
            warning(...
                'parse_keyword_args:unknown_key', ...
                'Unknown key in key-value list: %s, ignoring\n', key);
            continue
        elseif match_count > 1 ,
            warning(...
                'parse_keyword_args:repeated_key', ...
                'Repeated key in key-default list: %s, using first match\n', key);            
        end

        % If we get here, we know there's at least one match, so we take the first one
        matching_valid_key_index = matching_valid_key_indices(1) ;

        % Assign the value from this key-value pair to the right place in the output
        varargout{matching_valid_key_index} = value ;
    end
end
