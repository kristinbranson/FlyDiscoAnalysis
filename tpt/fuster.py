#!/usr/bin/env python

import time
import math
from tpt.utilities import *



def get_bjobs_lines(job_ids) :
    job_id_count = len(job_ids) 
    job_id_count_per_call = 10000 
    batch_count = math.ceil(job_id_count / job_id_count_per_call) 
    bjobs_lines_from_batch_index = [None] * batch_count
    for batch_index in range(batch_count) :
        first_job_index = job_id_count_per_call * batch_index
        last_job_index = min(job_id_count_per_call * (batch_index+1), job_id_count) 
        job_ids_this_batch = job_ids[first_job_index:last_job_index]
        job_ids_as_strings = [ ('%d' % job_id) for job_id in job_ids_this_batch ]
        command_line = ['bjobs'] + job_ids_as_strings
        (status, stdout) = run_subprocess_and_return_code_and_stdout(command_line) 
        if status != 0 :
            raise RuntimeError('There was a problem running the command "%s".  Return code: %d.  Stdout/stderr:\n%s' % (str(command_line), status, stdout) )
            
        lines = stdout.split('\n')   # Want a list of lines
        if len(lines)<1 :
            raise RuntimeError('There was a problem submitting the bjobs command "%s".  Unable to parse output.  Output was:\n%s', (str(command_line), stdout) ) 
        job_ids_this_batch_count = len(job_ids_this_batch) 
        bjobs_lines = lines[1:(job_ids_this_batch_count+1)]   # drop the header
        bjobs_lines_from_batch_index[batch_index] = bjobs_lines 
    result = flatten(bjobs_lines_from_batch_index)
    return result



def get_single_bsub_job_status(job_id) :
    # Possible results are {-1,0,+1}.
    #   -1 means errored out
    #    0 mean running or pending
    #   +1 means completed successfully
    
    if math.isfinite(job_id) :
        # This means the job has not been submitted yet
        result = math.nan 
    elif job_id == -1 :
        # This is a job that was run locally and exited cleanly
        result = +1 
    elif job_id == -2 :
        # This is a job that was run locally and errored
        result = -1 
    else :
        command_line = ['bjobs', str(job_id)]
        (status, stdout) = run_subprocess_and_return_code_and_stdout(command_line) 
        if status != 0 :
            raise RuntimeError('There was a problem running the command %s.  The return code was %d' % (command_line, status)) 
        lines = stdout.split('\n')
        if len(lines)<2 :
            raise RuntimeError('There was a problem submitting the bjobs command "%s".  Unable to parse output.  Output was:\n%s' % (space_out(command_line), stdout)) 
        bjobs_line = lines[1] 
        tokens = bjobs_line.split() 
        if len(tokens)<3 :
            raise RuntimeError('There was a problem submitting the bjobs command "%s".  Unable to parse output.  Output was:\n%s' % (space_out(command_line), stdout)) 
        lsf_status = tokens[2]   # Should be string like 'DONE', 'EXIT', 'RUN', 'PEND', etc.    
        if lsf_status == 'DONE' :
            result = +1 
        elif lsf_status=='EXIT' :
            # This seems to indicate an exit with something other than a 0 return code
            result = -1 
        elif lsf_status=='PEND' or lsf_status=='RUN' or lsf_status=='UNKWN' :
            result = 0 
        else :
            raise RuntimeError('Unknown bjobs status string: %s' % lsf_status) 
    return result



def get_bsub_job_status(job_ids) :
    # Possible results are {-1,0,+1,nan}.
    #   -1 means errored out
    #    0 mean running or pending
    #   +1 means completed successfully
    #  nan means the corresponding job_id was nan
    
    job_count = len(job_ids) 
    result = [math.nan] * job_count 
    has_not_been_submitted = list(map(math.isnan, job_ids))
    assign_where_true_bang(result, has_not_been_submitted, 0)
    was_run_locally = list(map(lambda job_id : (job_id<0), job_ids))   # means the job was run locally
    was_run_locally_and_exited_cleanly = list(map(lambda job_id : (job_id==-1), job_ids))
    was_run_locally_and_errored = list(map(lambda job_id : (job_id==-2), job_ids))
    assign_where_true_bang(result, was_run_locally_and_exited_cleanly, +1)
    assign_where_true_bang(result, was_run_locally_and_errored, -1)
    if all(elementwise_list_or(has_not_been_submitted, was_run_locally)) :
        return result
    was_submitted = [ (not el) for el in elementwise_list_or(was_run_locally, has_not_been_submitted) ]
    submitted_job_ids = ibb(job_ids, was_submitted) 
    bjobs_lines = get_bjobs_lines(submitted_job_ids) 
    bjobs_line_index = 0
    for job_index in range(job_count) :
        if not was_submitted[job_index] :
            continue 
        job_id = job_ids[job_index] 
        bjobs_line = bjobs_lines[job_index]
        tokens = bjobs_line.split()
        if len(tokens)<3 :
            raise RuntimeError('There was a problem with a bjobs command.  Unable to parse output.  Output was: %s' % bjobs_line) 
        running_job_id_as_string = tokens[0] 
        running_job_id = int(running_job_id_as_string) 
        if running_job_id != job_id :
            raise RuntimeError('The running job id (%d) doesn''t match the job id (%d)' % (running_job_id, job_id) ) 
        lsf_status = tokens[2]   # Should be string like 'DONE', 'EXIT', 'RUN', 'PEND', etc.
        if lsf_status == 'DONE' :
            running_job_status_code = +1 
        elif lsf_status == 'EXIT' :
            # This seems to indicate an exit with something other than a 0 return code
            running_job_status_code = -1 
        elif ( lsf_status=='PEND' or lsf_status=='RUN' or lsf_status=='UNKWN' or
               lsf_status=='SSUSP' or lsf_status=='PSUSP' or lsf_status=='USUSP' ) :
            running_job_status_code = 0 
        else :
            raise RuntimeError('Unknown bjobs status string: %s', lsf_status) 
        result[job_index] = running_job_status_code 
        bjobs_line_index = bjobs_line_index + 1 
    return result



def determine_which_jobs_to_submit(slot_count_from_submittable_index, maximum_slot_count) :
    # Determine which of the submittable jobs will be submitted, given how
    # many slots each submittable job needs, and the maximum number of slots we can use.
    submittable_count = len(slot_count_from_submittable_index) 
    will_submit_from_submittable_index = [False] * submittable_count
    slots_used_so_far = 0 
    for submittable_index in range(submittable_count) :
        slot_count_this_submittable = slot_count_from_submittable_index[submittable_index]
        putative_slots_used = slots_used_so_far + slot_count_this_submittable 
        if putative_slots_used <= maximum_slot_count :
            will_submit_from_submittable_index[submittable_index] = True 
            slots_used_so_far = putative_slots_used 
            if slots_used_so_far >= maximum_slot_count :
                break
    return will_submit_from_submittable_index



def bsub(command_line_as_list, do_actually_submit=True, slot_count=1, stdouterr_file_name='/dev/null', options_as_list=[]) :
    # Wrapper for LSF bsub command.  Returns job id as a double.
    # Throws error if anything goes wrong.
    if (stdouterr_file_name is None) or len(stdouterr_file_name)==0 :
        stdouterr_file_name = '/dev/null' 
    if do_actually_submit :
        bsub_command_line_as_string = \
            ( [ 'bsub', '-n', str(slot_count), '-oo', stdouterr_file_name, '-eo', stdouterr_file_name ] + 
              options_as_list + 
              command_line_as_list )      
        # printf('%s\n', bsub_command) 
        raw_stdout = run_subprocess_and_return_stdout(bsub_command_line_as_string)
        stdout = raw_stdout.strip()   # There are leading newlines and other nonsense in the raw version
        raw_tokens = stdout.split()
        is_token_nonempty = [ len(str)>0 for str in raw_tokens ]
        tokens = ibb(raw_tokens, is_token_nonempty) 
        is_job_token = [ token=='Job' for token in tokens ] 
        job_token_indices = where(is_job_token) 
        if isempty(job_token_indices) :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        job_token_index = job_token_indices[0]
        if len(tokens) < job_token_index+3 :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        if tokens[job_token_index+2] != 'is' :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        if tokens[job_token_index+3] != 'submitted' :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        job_id_token = tokens[job_token_index+1]
        if len(job_id_token)<2 :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        if job_id_token[0] != '<' or job_id_token[-1] != '>' :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
        job_id_as_string = job_id_token[1:-1]
        try :
            job_id = int(job_id_as_string) 
        except ValueError :
            raise RuntimeError('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s' %
                               (repr(bsub_command_line_as_string), stdout) )
    else :
        # Just call the def locally, but use a try/catch to make it more robust.
        try :
            run_subprocess_live(command_line_as_list, check=True)
            job_id = -1   # represents a job that was run locally and exited cleanly
        except RuntimeError as e :
            printf('Encountered an error while running a local bsub job.  Here''s some information about the error:\n') 
            printf('%s\n' % str(e)) 
            job_id = -2   # represents a job that was run locally and errored            
    return job_id



def update_job_status_from_job_index(old_job_status_from_job_index, job_id_from_job_index) :
    '''
    Calls get_bsub_job_status() on in-progress jobs to generate an updated job_status_from_job_index.
    Note that this is does *not* mutate the input at all.
    Ideally we could just call get_bsub_job_status() on all the jobs, but if the "conductor" jobs runs for a long time,
    bjobs will eventually throw errors because it has forgotten about long-completed jobs.
    '''
    was_in_progress_from_job_index = [ job_status==0 for job_status in old_job_status_from_job_index ]
    job_id_from_was_in_progress_job_index = ibb(job_id_from_job_index, was_in_progress_from_job_index)
    job_status_from_was_in_progress_job_index = get_bsub_job_status(job_id_from_was_in_progress_job_index)
    job_index_from_was_in_progress_job_index = where(was_in_progress_from_job_index)
    job_status_from_job_index = overlay_at(old_job_status_from_job_index, job_index_from_was_in_progress_job_index, job_status_from_was_in_progress_job_index)
    return job_status_from_job_index



class bqueue_type :    
    def __init__(self, do_actually_submit=True, maximum_running_slot_count=math.inf) :
        self._bsub_option_list_from_job_index = []
        self._command_line_as_list = [] 
        self._slot_count_from_job_index = []
        self._job_id_from_job_index = []
        self._has_been_submitted_from_job_index = []
        self._stdouterr_file_name_from_job_index = []
        self._do_actually_submit = do_actually_submit 
        self._maximum_running_slot_count = maximum_running_slot_count 
        self._job_status_from_job_index = []
        
    def queue_length(self) :
        result = len(self._has_been_submitted_from_job_index) 
        return result
    
    def enqueue(self, slot_count, stdouterr_file_name, bsub_options_as_list, command_line_as_list) :
        job_index = self.queue_length() + 1
        self._command_line_as_list.append(command_line_as_list)
        self._job_id_from_job_index.append(math.nan)
        self._has_been_submitted_from_job_index.append(False)            
        self._slot_count_from_job_index.append(slot_count)
        self._stdouterr_file_name_from_job_index.append(stdouterr_file_name)
        self._bsub_option_list_from_job_index.append(bsub_options_as_list)
        self._job_status_from_job_index.append(math.nan)
    
    def run(self, maximum_wait_time=math.inf, do_show_progress_bar=True) :
        # Possible job_statuses are {-1,0,+1,math.nan}.
        #   -1 means errored out
        #    0 mean running or pending
        #   +1 means completed successfully
        #   math.nan means not yet submitted
        
        have_all_exited = False 
        is_time_up = False 
        job_count = self.queue_length() 
        job_status_from_job_index = self._job_status_from_job_index
        if do_show_progress_bar :
            progress_bar = progress_bar_object(job_count) 
        ticId = tic() 
        while not have_all_exited and not is_time_up :
            old_job_status_from_job_index = job_status_from_job_index
            job_status_from_job_index = update_job_status_from_job_index(old_job_status_from_job_index, self._job_id_from_job_index)
            is_in_progress_from_job_index = [ job_status==0 for job_status in job_status_from_job_index ]
            carryover_slot_count = sum(ibb(self._slot_count_from_job_index, is_in_progress_from_job_index)) 
            maximum_new_slot_count = self._maximum_running_slot_count - carryover_slot_count 
            if maximum_new_slot_count > 0 :
                is_submittable_from_job_index = [ math.isnan(job_status) for job_status in job_status_from_job_index ]
                job_index_from_submittable_index = where(is_submittable_from_job_index) 
                slot_count_from_submittable_index = ibl(self._slot_count_from_job_index, job_index_from_submittable_index)
                will_submit_from_submittable_index = determine_which_jobs_to_submit(slot_count_from_submittable_index, maximum_new_slot_count) 
                job_indices_to_submit = ibb(job_index_from_submittable_index, will_submit_from_submittable_index) 
                jobs_to_submit_count = len(job_indices_to_submit) 
                for i in range(jobs_to_submit_count) :
                    job_index = job_indices_to_submit[i] 
                    command_line_as_list = self._command_line_as_list[job_index]
                    this_slot_count = self._slot_count_from_job_index[job_index]
                    this_stdouterr_file_name = self._stdouterr_file_name_from_job_index[job_index]
                    this_bsub_options_as_list = self._bsub_option_list_from_job_index[job_index]
                    this_job_id = \
                        bsub(command_line_as_list, 
                             self._do_actually_submit,
                             this_slot_count,
                             this_stdouterr_file_name,
                             this_bsub_options_as_list) 
                    self._job_id_from_job_index[job_index] = this_job_id 
                    if self._do_actually_submit :
                        job_status_from_job_index [job_index] = 0  # means running or pending 
                    else :
                        # This means the job was run locally, and has already either succeeded or failed
                        # In this case the value returned from bsub() indicates which.
                        if this_job_id == -1 :
                            job_status_from_job_index [job_index] = +1
                        elif this_job_id == -2 :
                            job_status_from_job_index [job_index] = -1
                        else :
                            raise RuntimeError('Programming error')
            has_job_exited = [ (job_status==-1 or job_status==+1) for job_status in job_status_from_job_index ]  
            exited_job_count = sum(has_job_exited)
            had_job_exited = [ (job_status==-1 or job_status==+1) for job_status in old_job_status_from_job_index ]  
            last_exited_job_count = sum(had_job_exited)
            newly_exited_job_count = exited_job_count - last_exited_job_count 
            if do_show_progress_bar :
                progress_bar.update(newly_exited_job_count) 
            have_all_exited = (exited_job_count==job_count) 
            self._job_status_from_job_index = job_status_from_job_index  # not necessary, but nice to keep things up to date
            if not have_all_exited :
                if self._do_actually_submit :
                    time.sleep(1) 
                is_time_up = (toc(ticId) > maximum_wait_time) 
        return job_status_from_job_index



def bwait(job_ids, maximum_wait_time=math.inf, do_show_progress_bar=True) :
    have_all_exited = False 
    is_time_up = False 
    job_count = len(job_ids) 
    if do_show_progress_bar :
        progress_bar = progress_bar_object(job_count) 
    last_exited_job_count = 0 
    ticId = tic() 
    while not have_all_exited and not is_time_up :
        job_statuses = get_bsub_job_status(job_ids) 
        has_job_exited = [status!=0 for status in job_statuses]
        exited_job_count = sum(has_job_exited) 
        newly_exited_job_count = exited_job_count - last_exited_job_count 
        if do_show_progress_bar :
            progress_bar.update(newly_exited_job_count) 
        have_all_exited = (exited_job_count==job_count)         
        if ~have_all_exited :
            time.sleep(10) 
            is_time_up = (toc(ticId) > maximum_wait_time) 
        last_exited_job_count = exited_job_count 
    return job_statuses



def test_bqueue() :
    do_actually_submit = True 
    max_running_slot_count = 5 
    bsub_options_as_list = [ '-P', 'scicompsoft', '-W', '59', '-J', 'test-bqueue' ]
    slots_per_job = 1 
    stdouterr_file_path = ''   # will go to /dev/null
    bqueue = bqueue_type(do_actually_submit, max_running_slot_count) 
    command_line_as_list = ['/usr/bin/sleep', '20']

    job_count = 10 
    for job_index in range(job_count) :
        bqueue.enqueue(slots_per_job, stdouterr_file_path, bsub_options_as_list, command_line_as_list)

    maximum_wait_time = 200 
    do_show_progress_bar = True 
    tic_id = tic() 
    job_statuses = bqueue.run(maximum_wait_time, do_show_progress_bar)
    elapsed_time = toc(tic_id)
    print('Elapsed time was %g seconds' % elapsed_time)
    print('Final job_statuses: %s' % str(job_statuses)) 
    if all([job_status==+1 for job_status in job_statuses]) :
        print('Test passed.')
    else:
        print('Test failed.')



# If called from command line, run the test(s)
if __name__ == "__main__":
    test_bqueue()
