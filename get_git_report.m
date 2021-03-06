function breadcrumb_string = get_git_report(source_repo_folder_path)
    original_pwd = pwd() ;
    cleaner = onCleanup(@()(cd(original_pwd))) ;
    cd(source_repo_folder_path) ;
    
    % Get the matalb version
    matlab_ver_string = version() ;
    
    % This is hard to get working in a way that overrides
    % 'url."git@github.com:".insteadOf https://github.com/' for a single command.
    % Plus it hits github every time you run, which seems fragile...
    % % Make sure the git remote is up-to-date
    % system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git remote update') ;    
    
    % Get the git hash
    stdout = system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git rev-parse --verify HEAD') ;
    commit_hash = strtrim(stdout) ;

    % Get the git remote report
    git_remote_report = system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git remote -v') ;    
    
    % Get the git status
    git_status = system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git status') ;    
    
    % Get the recent git log
    git_log = system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git log --graph --oneline --max-count 10 | cat') ;
        
    % Write a file with the commit hash into the folder, for good measure
    breadcrumb_string = sprintf('Matlab version:\n%s\n\nSource repo:\n%s\n\nCommit hash:\n%s\n\nRemote info:\n%s\n\nGit status:\n%s\n\nGit log:\n%s\n\n', ...
                                matlab_ver_string, ...
                                source_repo_folder_path, ...
                                commit_hash, ...
                                git_remote_report, ...
                                git_status, ...
                                git_log) ;
end
