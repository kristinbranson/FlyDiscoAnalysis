#!/usr/bin/python3

import sys
import os
import subprocess



class cd:
    """Context manager for changing the current working directory, and automagically changing back when done"""
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.old_path)



def run_subprocess_live(command_as_list, check=True, shell=False) :
    '''
    Call an external executable, with live display of the output.
    '''
    p = subprocess.Popen(command_as_list, shell=shell)
    return_code = p.wait()
    if check :
        if return_code != 0 :
            raise RuntimeError("Running %s returned a non-zero return code: %d" % (str(command_as_list), return_code))
    return return_code



def run_subprocess_and_return_stdout(command_as_list, shell=False) :
    completed_process = \
        subprocess.run(command_as_list, 
                       stdout=subprocess.PIPE,
                       encoding='utf-8',
                       check=True, 
                       shell=shell)
    stdout = completed_process.stdout
    #print('Result: %s' % result)                   
    return stdout



def get_git_report(source_repo_folder_path) :
    with cd(source_repo_folder_path) as _ :
        # Get the Python version
        python_ver_string = sys.version
        
        # This is hard to get working in a way that overrides
        # 'url."git@github.com:".insteadOf https://github.com/' for a single command.
        # Plus it hits github every time you run, which seems fragile\
        # % Make sure the git remote is up-to-date
        # system_with_error_handling('env GIT_SSL_NO_VERIFY=true GIT_TERMINAL_PROMPT=0 git remote update')     
        
        # Get the git hash
        so = run_subprocess_and_return_stdout(
                ['/usr/bin/env', 'GIT_SSL_NO_VERIFY=true', 'GIT_TERMINAL_PROMPT=0', '/usr/bin/git', 'rev-parse', '--verify', 'HEAD'])
        commit_hash = so.strip()

        # Get the git remote report
        git_remote_report = run_subprocess_and_return_stdout(
            ['/usr/bin/env', 'GIT_SSL_NO_VERIFY=true', 'GIT_TERMINAL_PROMPT=0', '/usr/bin/git',  'remote',  '-v'])     

        # Get the git status
        git_status = run_subprocess_and_return_stdout(
            ['/usr/bin/env', 'GIT_SSL_NO_VERIFY=true', 'GIT_TERMINAL_PROMPT=0', '/usr/bin/git', 'status'])     
        
        # Get the recent git log
        git_log = run_subprocess_and_return_stdout(
            ['/usr/bin/env', 'GIT_SSL_NO_VERIFY=true', 'GIT_TERMINAL_PROMPT=0', '/usr/bin/git', 'log', '--graph', '--oneline', '--max-count', '10']) 
            
    # Package everything up into a string
    breadcrumb_string = 'Python version:\n%s\n\nSource repo:\n%s\n\nCommit hash:\n%s\n\nRemote info:\n%s\n\nGit status:\n%s\n\nGit log:\n%s\n\n' % \
                        (python_ver_string, 
                         source_repo_folder_path, 
                         commit_hash, 
                         git_remote_report, 
                         git_status, 
                         git_log) 

    return breadcrumb_string



def transfero_FlyDiscoPipeline_wrapper_wrapper(raw_experiment_folder_path):
    experiment_folder_path = os.path.realpath(os.path.abspath(raw_experiment_folder_path))

    # Find the matlab code
    this_script_path = os.path.realpath(__file__)
    this_script_folder_path = os.path.dirname(this_script_path)

    ## List all the environment variables
    #print('Environment:')
    #for name, value in os.environ.items():
    #    print("%s: %s" % (name, value))
    #print('')
        
    # Output the git report for this folder
    git_report = get_git_report(this_script_folder_path) 
    print(git_report) 
    
    # Flush stdout output now, so it comes before an stdout output of the subprocess
    sys.stdout.flush()

    with cd(this_script_folder_path) as _ :
        matlab_command_line = \
            "modpath; options = cell(1, 0); transfero_FlyDiscoPipeline_wrapper('%s', options)" % \
                experiment_folder_path
        #print("Matlab command line is: %s" % matlab_command_line)
        command_line_as_list = ['/usr/bin/xvfb-run', '-d', '/misc/local/matlab-2019a/bin/matlab', '-batch', matlab_command_line]
        #print("Subprocess command line as list is: %s" % repr(command_line_as_list))        
        sys.stdout.flush()  # Flush stdout output now, so it comes before an stdout output of the subprocess
        run_subprocess_live(command_line_as_list)



if __name__ == "__main__":
    if len(sys.argv)==2 :
        transfero_FlyDiscoPipeline_wrapper_wrapper(sys.argv[1])
    else:
        raise RuntimeError('Need exactly one argument')
