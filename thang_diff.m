function thang_diff(test, ref)
    % Recursively enumerate difference between structures, cell arrays, whatever.
    thang_diff_helper(test, ref, 'test', 'ref') ;
end



function thang_diff_helper(test, ref, test_prefix, ref_prefix)
    if isequaln(test, ref) ,
        return
    end    
    if isnumeric(ref) ,
        if ~isequal(class(test), class(ref)) ,
            fprintf('%s is of class %s, %s is of class %s\n', ref_prefix, class(ref), test_prefix, class(test)) ;
            return
        end
        % both numeric, of same class
        ref_size = size(ref) ;
        test_size = size(test) ;
        if ~isequaln(test_size, ref_size) ,
            fprintf('%s is a %s array of size %s, %s the same but of size %s\n', ref_prefix, class(ref), disp(ref_size), test_prefix, disp(test_size)) ;
            return
        end
        % both numeric, of same class, same size 
        abs_diff = abs(test-ref) ;
        max_abs_diff = max(abs_diff, [], 'all') ;
        min_abs_diff = min(abs_diff, [], 'all') ;
        fprintf('%s and %s are both numeric, but are not the same.  max abs diff = %g, min abs diff = %g\n', ...
                ref_prefix, test_prefix, max_abs_diff, min_abs_diff) ;
        return
    elseif isstruct(ref) ,
        if ~isstruct(test) ,
            fprintf('%s is a struct, %s is not\n', ref_prefix, test_prefix) ;
            return
        end
        % test is a struct
        ref_size = size(ref) ;
        test_size = size(test) ;
        if ~isequaln(test_size, ref_size) ,
            fprintf('%s is a struct of size %s, %s is a struct of size %s\n', ref_prefix, disp(ref_size), test_prefix, disp(test_size)) ;
            return
        end
        ref_field_names = fieldnames(ref) ;
        test_field_names = fieldnames(test) ;
        extra_field_names = setdiff(ref_field_names, test_field_names) ;
        missing_field_names = setdiff(test_field_names, ref_field_names) ;
        if ~isempty(extra_field_names) || ~isempty(missing_field_names) ,
            fprintf('%s has different field names than %s\n', test_prefix, ref_prefix) ;
            return
        end        
        ref_count = numel(ref) ;
        if ref_count == 0 ,
            return 
        elseif ref_count==1 ,
            scalar_struct_diff(test, ref, test_prefix, ref_prefix)
        else
            for i = 1 : ref_count ,
                sub_ref = ref(i) ;
                sub_test = test(i) ;
                sub_ref_prefix = sprintf('%s(%d)', ref_prefix, i) ;
                sub_test_prefix = sprintf('%s(%d)', test_prefix, i) ;
                scalar_struct_diff(sub_test, sub_ref, sub_test_prefix, sub_ref_prefix)
            end
        end
    elseif iscell(ref)
        if ~iscell(test) ,
            fprintf('%s is a cell array, %s is not\n', ref_prefix, test_prefix) ;
            return
        end
        % test is a cell array
        ref_size = size(ref) ;
        test_size = size(test) ;
        if ~isequaln(test_size, ref_size) ,
            fprintf('%s is a cell array of size %s, %s is a cell array of size %s\n', ref_prefix, disp(ref_size), test_prefix, disp(test_size)) ;
            return
        end
        ref_count = numel(ref) ;
        if ref_count == 0 ,
            return 
        else
            for i = 1 : ref_count ,
                sub_ref = ref{i} ;
                sub_test = test{i} ;
                sub_ref_prefix = sprintf('%s{%d}', ref_prefix, i) ;
                sub_test_prefix = sprintf('%s{%d}', test_prefix, i) ;
                thang_diff_helper(sub_test, sub_ref, sub_test_prefix, sub_ref_prefix)
            end
        end        
    elseif ischar(ref)
        if ~ischar(test) ,
            fprintf('%s is a char array, %s is not\n', ref_prefix, test_prefix) ;
            return
        end
        fprintf('%s and %s are both char arrays, but are not the same\n', ref_prefix, test_prefix) ;
        return        
    else
        error('%s is of class %s, which is not supported', ref_prefix, class(ref)) ;
    end
end



function scalar_struct_diff(s_test, s_ref, test_prefix, ref_prefix)
    % s_test, s_ref should be scalar structs with the same field names
    if isequaln(s_test, s_ref) ,
        return
    end
    ref_field_names = fieldnames(s_ref) ;
    field_count = length(ref_field_names) ;
    for i = 1 : field_count ,
        field_name = ref_field_names{i} ;
        % use getfield() b/c it works even in the case of field names like '1'
        sub_ref = getfield(s_ref, field_name) ;  %#ok<GFLD>
        sub_test = getfield(s_test, field_name) ;  %#ok<GFLD>
        sub_ref_prefix = sprintf('%s.%s', ref_prefix, field_name) ;
        sub_test_prefix = sprintf('%s.%s', test_prefix, field_name) ;
        thang_diff_helper(sub_test, sub_ref, sub_test_prefix, sub_ref_prefix) ;
    end
end
