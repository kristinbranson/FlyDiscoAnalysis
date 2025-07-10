import sys
import os
import pwd
import time
import subprocess
import io
import datetime
import ast
import copy
import shlex
import stat
import tempfile



class cd:
    """Context manager for changing the current working directory, and automagically changing back when done"""
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.old_path)



def boolean_from_string(s) :
    if s.lower()=='true' :
        return True
    elif s.lower()=='false' :
        return False
    else :
        raise RuntimeError('Unable to convert string "%s" to True or False')



def tic() :
    return time.time()



def toc(t) :
    return time.time() - t



def listmap(f, lst) :
    return list(map(f, lst))



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



def printf(*args, **kwargs) :
    '''
    print() but without the newline.
    '''
    print(*args, end='', **kwargs)



def printe(*args, **kwargs):
    '''
    print() but to stderr.
    '''
    print(*args, file=sys.stderr, **kwargs)



def printfe(*args, **kwargs):
    '''
    print() but to stderr, and without the newline.
    '''
    print(*args, file=sys.stderr, end='', **kwargs)



def run_subprocess_and_return_stdout(command_as_list, shell=False) :
    completed_process = \
        subprocess.run(command_as_list, 
                       stdout=subprocess.PIPE,
                       encoding='utf-8',
                       check=False, 
                       shell=shell)
    stdout = completed_process.stdout    
    return_code = completed_process.returncode
    if return_code != 0 :
        raise RuntimeError('Command %s returned nonzero return code %d.\nstdout:\n%s\n' 
                           % (str(command_as_list), return_code, stdout) )
    return stdout



def run_subprocess_and_return_stdout_and_stderr(command_as_list, shell=False) :
    completed_process = \
        subprocess.run(command_as_list, 
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, 
                       encoding='utf-8',
                       check=False, 
                       shell=shell)
    stdout = completed_process.stdout    
    stderr = completed_process.stderr    
    return_code = completed_process.returncode
    if return_code != 0 :
        raise RuntimeError('Command %s returned nonzero return code %d.\nstdout:\n%s\nstderr:\n%s\n' 
                           % (str(command_as_list), return_code, stdout, stderr) )
    return (stdout, stderr)



def run_subprocess_and_return_code_and_stdout(command_as_list, shell=False) :
    completed_process = \
        subprocess.run(command_as_list, 
                       stdout=subprocess.PIPE,
                       encoding='utf-8',
                       check=False, 
                       shell=shell)
    stdout = completed_process.stdout
    return_code = completed_process.returncode
    #print('Result: %s' % result)                   
    return (return_code, stdout)



def run_subprocess_and_return_code_and_stdout_and_stderr(command_as_list, shell=False) :
    completed_process = \
        subprocess.run(command_as_list, 
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       encoding='utf-8',
                       check=False, 
                       shell=shell)
    stdout = completed_process.stdout
    stderr = completed_process.stderr    
    return_code = completed_process.returncode
    #print('Result: %s' % result)                   
    return (return_code, stdout, stderr)



def run_subprocess_and_return_code(command_as_list, shell=False) :
    '''
    Run the subprocess, with stdout/stderr going to those of the parent process.
    Return the return code.  *Don't* throw an exception for a nonzero return code.
    '''
    completed_process = \
        subprocess.run(command_as_list, 
                       check=False, 
                       shell=shell)
    return_code = completed_process.returncode
    #print('Result: %s' % result)                   
    return return_code



def run_subprocess(command_as_list, shell=False) :
    '''
    Run the subprocess, with stdout/stderr going to those of the parent process.
    Throw an exception if there's a problem, otherwise return without 
    returning anything.
    '''
    completed_process = \
        subprocess.run(command_as_list, 
                       check=True, 
                       shell=shell)



def run_subprocess_live_and_return_stdouterr(command_as_list, check=True, shell=False) :
    '''
    Call an external executable, with live display of the output.  
    Return stdout+stderr as a string.
    '''
    with subprocess.Popen(command_as_list, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.STDOUT, 
                          bufsize=1, 
                          encoding='utf-8', 
                          shell=shell) as p, io.StringIO() as buf:
        for line in p.stdout :
            print(line, end='')
            buf.write(line)
        p.communicate()  # Seemingly needed to make sure returncode is set.  
                         # Hopefully will not deadlock b/c we've already 
                         # exhausted stdout.
        return_code = p.returncode
        if check :
            if return_code != 0 :
                raise RuntimeError("Running %s returned a non-zero return code: %d" % (str(command_as_list), return_code))
        stdouterr = buf.getvalue()     
    return (stdouterr, return_code)



def run_subprocess_live(command_as_list, check=True, shell=False) :
    '''
    Call an external executable, with live display of the output.
    '''
    with subprocess.Popen(command_as_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, encoding='utf-8', shell=shell) as p, io.StringIO() as buf:
        for line in p.stdout :
            print(line, end='')
            buf.write(line)
        p.communicate()  # Seemingly needed to make sure returncode is set.  Hopefully will not deadlock b/c we've already exhausted stdout.    
        return_code = p.returncode
        if check :
            if return_code != 0 :
                raise RuntimeError("Running %s returned a non-zero return code: %d" % (str(command_as_list), return_code))
    return return_code



def run_subprocess_with_log_and_return_code(command_as_list, log_file_name, shell=False) :
    '''
    Call an external executable, with stdout+stderr to log file.
    '''
    with open(log_file_name, 'w') as fid:
        completed_process = subprocess.run(command_as_list, stdout=fid, stderr=subprocess.STDOUT, encoding='utf-8', shell=shell, check=False)
        return_code = completed_process.returncode
    return return_code



def space_out(lst) :
    '''
    Given a list of strings, return a single string with the list concatenated, but with spaces between them.
    '''
    result = '' 
    count = len(lst) 
    for i in range(count) :
        if i==0 :
            result = lst[i] 
        else :
            result = result + ' ' + lst[i] 
    return result



def get_user_name():
    return pwd.getpwuid(os.getuid())[0]



def listsetdiff(lst1, lst2) :
    # Set difference for lists
    s2 = set(lst2)
    result = [el for el in lst1 if el not in s2]    
    return result



def listsetintersect(lst1, lst2) :
    # Set intersection for lists
    s2 = set(lst2)
    result = [el for el in lst1 if el in s2]    
    return result



def ibl(lst, old_index_from_new_index) :
    # "Index By List"
    # Designed to mimic Matlab's x(old_index_from_new_index) syntax when old_index_from_new_index is an array of indices
    new_count = len(old_index_from_new_index)
    old_count = len(lst)
    if new_count>old_count :
        raise RuntimeError('old_index_from_new_index has more elements (%d) than lst (%d)' % (new_count, old_count))
    result = [None] * new_count  # pre-allocate output, although I'm not sure this helps with perf
    for new_index in range(new_count) :
        old_index = old_index_from_new_index[new_index]
        result[new_index] = lst[old_index]
    return result



def ibb(lst, pred) :
    # "Index By Boolean"
    # Designed to mimic Matlab's x(is_something) syntax when is_something is a boolan array
    result =[]
    for i in range(len(lst)) :
        if pred[i]:
            result.append(lst[i])
    return result



def ibbn(lst, pred) :
    # "Index By Boolean Negated"
    # Designed to mimic Matlab's x(~is_something) syntax when is_something is a boolan array
    result =[]
    for i in range(len(lst)) :
        if not pred[i]:
            result.append(lst[i])
    return result



def list_fif(pred, if_true, if_false) :
    # If pred[i] is true, result[i] is set to if_true[i], otherwise result[i] is set to if_false[i].
    n = len(pred)
    result = [None] * n
    for i in range(n) :
        if pred[i]:
            result[i] = if_true[i]
        else :
            result[i] = if_false[i]
    return result



def assign_where_true_bang(lst, pred, el) :
    # If pred[i] is true, lst[i] is set to el.
    # This mutates lst.
    for i in range(len(lst)) :
        if pred[i]:
            lst[i] = el



def overlay_at(lst, index_from_other_index, new_value_from_other_index) :
    '''
    Replaces the values in lst at indices given by index_from_other_index with values taken from new_value_from_other_index.
    lst is not mutated.
    '''
    other_index_count = len(new_value_from_other_index)
    result = copy.deepcopy(lst)
    for other_index in range(other_index_count) :
        index = index_from_other_index[other_index]
        result[index] = new_value_from_other_index[other_index]
    return result



def elementwise_list_and(a, b) :
    return [ (el_a and el_b) for (el_a,el_b) in zip(a, b) ]



def elementwise_list_or(a, b) :
    return [ (el_a or el_b) for (el_a,el_b) in zip(a, b) ]



def elementwise_list_not(a) :
  return [ (not el_a) for el_a in a ]



def flatten(lst):
    return [item for sublist in lst for item in sublist]



def argfilter(pred, lst) :
    # Like filter(), but returns the indices of elements that would be returned by filter(), not the elements themselves
    result = []
    n = len(lst)
    for i in range(n) :
        el = lst[i]
        if pred(el) :
            result.append(i)
    return result



def isempty(lst) :
    return (len(lst)==0)



def isladen(lst) :
    return (len(lst)>0)



class spinner_object :
    # Class for indicting progress
    def __init__(self, type="") :
        self.cursors = '|/-\\' 
        self.cursor_count = len(self.cursors) 
        self.cursor_index = 0 
        self.is_first_call = True 
        if type == 'mute' :
            self.is_mute = True 
        else :
            self.is_mute = False 

    def spin(self) :
        if not self.is_mute :
            if self.is_first_call :
                cursor = self.cursors[self.cursor_index]
                printfe(cursor) 
                self.is_first_call = False 
            else :
                printfe('\b')
                self.cursor_index = (self.cursor_index + 1) 
                if self.cursor_index >= self.cursor_count :
                    self.cursor_index = 0 
                cursor = self.cursors[self.cursor_index]
                printfe(cursor)

    def print(self, *varargin) :
        # Want things printed during spinning to look nice
        if not self.is_mute :
            printfe('\b\n')   # Delete cursor, then newline
            printfe(*varargin)   # print whatever
            cursor = self.cursors[self.cursor_index]   # get the same cursor back
            printfe(cursor)   # write it again on its own line

    def stop(self) :
        if not self.is_mute :
            printfe('\bdone.\n')



class progress_bar_object :
    # properties
    #     n_
    #     i_
    #     percent_as_displayed_last_
    #     did_print_at_least_one_line_
    #     did_print_final_newline_
    #     data_queue_
    # end
    
    def __init__(self, n) :
        self.n_ = n 
        self.i_ = 0 
        self.percent_as_displayed_last_ = [] 
        self.did_print_at_least_one_line_ = False 
        self.did_print_final_newline_ = False 
    
    def update(self, di=1) :
        # Should be called from within the for loop
        if self.did_print_final_newline_ :
            return
        self.i_ = min(self.i_ + di, self.n_) 
        i = self.i_ 
        n = self.n_ 
        percent = 100 if n==0 else 100*(i/n) 
        percent_as_displayed = round(percent*10)/10
        if percent_as_displayed != self.percent_as_displayed_last_ :
            if self.did_print_at_least_one_line_ :
                delete_bar = ''.join(['\b'] * (1+50+1+2+4+1)) 
                printfe(delete_bar) 
            bar = ''.join(['*'] * round(percent/2))
            printfe('[%-50s]: %4.1f%%' % (bar, percent_as_displayed)) 
            self.did_print_at_least_one_line_ = True 
        if i==n :
            if not self.did_print_final_newline_ :
                printfe('\n') 
                self.did_print_final_newline_ = True 
        self.percent_as_displayed_last_ = percent_as_displayed 




def where(t) :
    return [i for i, x in enumerate(t) if x]



def noop() :
    pass



def aware_datetime_from_timestamp(timestamp) :
    naive_datetime = datetime.datetime.fromtimestamp(timestamp)
    return naive_datetime.astimezone()  # aware timezone, represented in local TZ



def simple_dir(folder_name) :
    name_from_index = os.listdir(folder_name)
    path_from_index = list(map(lambda name: os.path.join(folder_name, name), 
                               name_from_index))
    is_folder_from_index = list(map(lambda path: os.path.isdir(path), 
                                    path_from_index))
    byte_count_from_index = list(map(lambda path: os.path.getsize(path),
                                     path_from_index))
    timestamp_from_index = list(map(lambda path: os.path.getmtime(path),
                                    path_from_index))
    datetime_from_index = list(map(aware_datetime_from_timestamp, timestamp_from_index))
    return (name_from_index, is_folder_from_index, byte_count_from_index, datetime_from_index)



def read_yaml_file_badly(file_name) :
    result = {}
    with open(file_name, 'r', encoding='UTF-8') as file:
        line_from_line_index = file.readlines()
    line_count = len(line_from_line_index)
    for line_index in range(line_count) :
        line = line_from_line_index[line_index].strip()
        index_of_colon = line.find(':')
        if index_of_colon < 0 :
            continue
        key = line[:index_of_colon].strip()
        value_as_string = line[index_of_colon+1:].strip()
        # print("value_as_string: %s" % value_as_string)
        value = ast.literal_eval(value_as_string)
        result[key] = value
    return result
    


class LockFile:
    '''
    Simple lock file implementation.  Definitely has the potential for race conditions, but
    only anticipate Transfero being run every 24 hours, so should be ok.
    '''
    def __init__(self, file_name):
        self._file_name = file_name 
        self._have_lock = False

    def __enter__(self):
        file_name = self._file_name
        if os.path.exists(file_name) :
            pass
        else :
            open(file_name, 'w').close()
            self._have_lock = True
        return self

    def __exit__(self, type, value, tb):
        if self._have_lock :
            file_name = self._file_name
            os.remove(file_name)

    def have_lock(self):
        return self._have_lock    



def run_remote_subprocess_and_return_stdout(user_name, host_name, remote_command_line_as_list) :
    '''
    Run the system command, but taking a list of tokens rather than a string, and
    running on a remote host.  Uses ssh, which needs to be set up for passowrdless
    login as the indicated user.
    Each element of command_line_as_list is escaped for bash, then composed into a
    single string, then submitted to system_with_error_handling().
    '''

    # Escape all the elements of command_line_as_list
    escaped_remote_command_line_as_list = [shlex.quote(el) for el in remote_command_line_as_list] 
    
    # Build up the command line by adding space between elements
    remote_command_line = space_out(escaped_remote_command_line_as_list)

    # Command line
    command_line_as_list = ['ssh', '-l', user_name, host_name, remote_command_line] ; 
    
    # Actually run the command
    stdout = run_subprocess_and_return_stdout(command_line_as_list)
    return stdout



def error_if_uncommited_changes(repo_path) :
    with cd(repo_path) as _ :
        stdout = run_subprocess_and_return_stdout(['git', 'status', '--porcelain=v1']) 
        trimmed_stdout = stdout.strip()  # Will be empty if no uncomitted changes
        if isladen(trimmed_stdout) :
            raise RuntimeError('The git repo seems to have uncommitted changes:\n%s' % stdout) 



def copy_local_repository_to_single_user_account(user_name, repository_folder_path):
    # Copy the folder over
    host_name = 'login2'   # Why not?
    repository_name = os.path.basename(repository_folder_path)
    printf('Copying %s into the %s user account...' % (repository_name, user_name) )
    run_remote_subprocess_and_return_stdout(user_name, host_name, ['rm', '-rf', repository_name]) 
    run_remote_subprocess_and_return_stdout(user_name, host_name, ['cp', '-R', '-T', repository_folder_path, repository_name]) 
    printf('done.\n') 



def clone_and_copy_github_repository_into_user_home_folders(url, username_from_user_index, branch_name='main') :
    # Get the repo name
    repository_name = os.path.basename(url)

    # Determine the folder path of this script
    this_script_path = os.path.realpath(__file__)
    script_folder_path = os.path.dirname(this_script_path)

    # Create a temporary folder
    with tempfile.TemporaryDirectory(dir=script_folder_path) as temp_folder_path :
        # Make the temp folder world-readable
        os.chmod(temp_folder_path, 
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

        with cd(temp_folder_path) as _ :
            # Clone the repo
            run_subprocess_and_return_code_and_stdout(
                ['git', 'clone', '-b', branch_name, '--recurse-submodules', url] )

            # Determine the cloned repo folder path
            repository_folder_path = os.path.join(temp_folder_path, repository_name)
    
            # Copy into all the given user account home folders
            for user_index in range(len(username_from_user_index)) :
                username = username_from_user_index[user_index]
                copy_local_repository_to_single_user_account(username, repository_folder_path)

            # If get here, everything went well
            printf('Successfully copied %s into all the *lab/*robot user accounts\n' % repository_name) 
