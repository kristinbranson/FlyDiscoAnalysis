function breadcrumb_string = get_git_report(source_repo_folder_path)
    original_pwd = pwd() ;
    cleaner = onCleanup(@()(cd(original_pwd))) ;
    cd(source_repo_folder_path) ;
    
    % Get the matalb version
    matlab_ver_string = version() ;
    
    % Make sure the git remote is up-to-date
    system_with_error_handling('git remote update') ;    
    
    % Get the git hash
    stdout = system_with_error_handling('git rev-parse --verify HEAD') ;
    commit_hash = strtrim(stdout) ;

    % Get the git remote report
    git_remote_report = system_with_error_handling('git remote -v') ;    
    
    % Get the git status
    git_status = system_with_error_handling('git status') ;    
    
    % Get the recent git log
    git_log = system_with_error_handling('git log --graph --oneline --max-count 10 | cat') ;
        
    % Write a file with the commit hash into the folder, for good measure
    breadcrumb_string = sprintf('Matlab version:\n%s\n\nSource repo:\n%s\n\nCommit hash:\n%s\n\nRemote info:\n%s\n\nGit status:\n%s\n\nGit log:\n%s\n\n', ...
                                matlab_ver_string, ...
                                source_repo_folder_path, ...
                                commit_hash, ...
                                git_remote_report, ...
                                git_status, ...
                                git_log) ;
end
