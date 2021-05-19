function error_if_uncommited_changes(repo_path)
    original_pwd =pwd() ;
    cleaner = onCleanup(@()(cd(original_pwd))) ;
    cd(repo_path) ;
    stdout = system_with_error_handling('git status --porcelain=v1') ;
    trimmed_stdout = strtrim(stdout) ;  % Will be empty if no uncomitted changes
    if ~isempty(trimmed_stdout) ,
        error('The git repo seems to have uncommitted changes:\n%s', stdout) ;
    end
    %stdout = system_with_error_handling('git rev-parse --verify HEAD') ;
    %result = strtrim(stdout) ;
end
