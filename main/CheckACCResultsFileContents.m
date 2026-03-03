function do_continue_pipeline = CheckACCResultsFileContents(expdir, dataloc_params)
    acc_results_as_text_file_name = dataloc_params.automaticcheckscompleteresultsfilestr ;
    acc_results_as_text_file_path = fullfile(expdir, acc_results_as_text_file_name) ;
    if ~exist(acc_results_as_text_file_path, 'file')
        fprintf('ACC results file %s does not exist\n', acc_results_as_text_file_name) ;
        do_continue_pipeline = false ;
        return
    end
    acc_results = ReadParams(acc_results_as_text_file_path) ;
    if isfield(acc_results, 'automated_pf') ,
        if strcmp(acc_results.automated_pf, 'P'),
            do_continue_pipeline = true ;
        elseif strcmp(acc_results.automated_pf, 'F') ,
            fprintf('ACC results file %s has automated_pf==F\n', acc_results_as_text_file_name) ;
            do_continue_pipeline = false ;
        else
            fprintf('ACC results file %s has automated_pf value that is not P nor F\n', acc_results_as_text_file_name) ;
            do_continue_pipeline = false ;
        end
    else
        fprintf('ACC results file %s has no automated_pf field\n', acc_results_as_text_file_name) ;
        do_continue_pipeline = false ;
    end
end
