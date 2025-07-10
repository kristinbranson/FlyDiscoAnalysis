classdef bqueue_type < handle
    properties (SetAccess = private)
        do_actually_submit = true
        do_try = true 
        bsub_option_string_from_job_index = cell(1,0) ;
        function_handle = cell(1,0)
        other_arguments = cell(1,0)
        %job_status = nan(1,0) ;
        maximum_running_slot_count = inf        
        slot_count_from_job_index = zeros(1,0) ;
        job_ids = zeros(1,0)
        has_job_been_submitted = false(1,0)
        %in_use_slot_count = 0 ;
        stdouterr_file_name_from_job_index = cell(1,0) ;
        do_use_xvfb = false ;
        submit_host_name = '' ;
    end
    
    methods
        function self = bqueue_type(do_actually_submit, do_try, maximum_running_slot_count, do_use_xvfb, submit_host_name)
            if ~exist('do_actually_submit', 'var') || isempty(do_actually_submit) ,
                do_actually_submit = true ;
            end
            if ~exist('do_try', 'var') || isempty(do_try) ,
                do_try = true ;
            end
            if ~exist('maximum_running_slot_count', 'var') || isempty(maximum_running_slot_count) ,
                maximum_running_slot_count = 400 ;
            end
            if ~exist('do_use_xvfb', 'var') || isempty(do_use_xvfb) ,
                do_use_xvfb = false ;
            end
            if ~exist('submit_host_name', 'var') || isempty(submit_host_name) ,                
                submit_host_name = '' ;
            end
            self.do_actually_submit = do_actually_submit ;
            self.do_try = do_try ;
            self.maximum_running_slot_count = maximum_running_slot_count ;
            self.do_use_xvfb = do_use_xvfb ;
            self.submit_host_name = submit_host_name ;
        end
        
        function result = queue_length(self)
            result = length(self.has_job_been_submitted) ;
        end
        
        function enqueue(self, slot_count, stdouterr_file_name, bsub_option_string, function_handle, varargin)
            job_index = self.queue_length() + 1 ;
            self.function_handle{1,job_index} = function_handle ;
            self.other_arguments(1,job_index) = {varargin(:)} ;            
            %self.job_status(1,job_index) = nan ;  % nan means unsubmitted
            self.job_ids(1,job_index) = nan ;            
            self.has_job_been_submitted(1,job_index) = false ;            
            self.slot_count_from_job_index(job_index) = slot_count ;
            self.stdouterr_file_name_from_job_index{1,job_index} = stdouterr_file_name ;
            self.bsub_option_string_from_job_index{1,job_index} = bsub_option_string ;
        end
        
        function job_statuses = run(self, maximum_wait_time, do_show_progress_bar)
            % Possible job_statuses are {-1,0,+1}.
            %   -1 means errored out
            %    0 mean running or pending
            %   +1 means completed successfully            
            if ~exist('maximum_wait_time', 'var') || isempty(maximum_wait_time) ,
                maximum_wait_time = inf ;
            end
            if ~exist('do_show_progress_bar', 'var') || isempty(do_show_progress_bar) ,
                do_show_progress_bar = true ;
            end

            have_all_exited = false ;
            is_time_up = false ;
            job_count = self.queue_length() ;
            if do_show_progress_bar ,
                progress_bar = progress_bar_object(job_count) ;
            end
            last_exited_job_count = 0 ;
            ticId = tic() ;
            while ~have_all_exited && ~is_time_up ,
                old_job_ids = self.job_ids ;                
                job_statuses = get_bsub_job_status(old_job_ids, self.submit_host_name) ; 
                is_job_in_progress = self.has_job_been_submitted & (job_statuses==0) ;
                carryover_slot_count = sum(self.slot_count_from_job_index(is_job_in_progress)) ;
                maximum_new_slot_count = self.maximum_running_slot_count - carryover_slot_count ;
                if maximum_new_slot_count > 0 ,
                    is_submittable_from_job_index = ~self.has_job_been_submitted ;
                    job_index_from_submittable_index = find(is_submittable_from_job_index) ;
                    slot_count_from_submittable_index = self.slot_count_from_job_index(job_index_from_submittable_index) ;
                    will_submit_from_submittable_index = determine_which_jobs_to_submit(slot_count_from_submittable_index, maximum_new_slot_count) ;
                    job_indices_to_submit = job_index_from_submittable_index(will_submit_from_submittable_index) ;
                    jobs_to_submit_count = length(job_indices_to_submit) ;
                    for i = 1 : jobs_to_submit_count ,
                        job_index = job_indices_to_submit(i) ;
                        this_function_handle = self.function_handle{job_index} ;
                        this_other_arguments = self.other_arguments{job_index} ;
                        this_slot_count = self.slot_count_from_job_index(job_index) ;
                        this_stdouterr_file_name = self.stdouterr_file_name_from_job_index{job_index} ;
                        this_bsub_option_string = self.bsub_option_string_from_job_index{job_index} ;
                        this_job_id = ...
                            bsub(self.do_actually_submit, ...
                                 self.do_try, ...
                                 this_slot_count, ...
                                 this_stdouterr_file_name, ...
                                 this_bsub_option_string, ...
                                 self.do_use_xvfb, ...
                                 self.submit_host_name, ...
                                 this_function_handle, ...
                                 this_other_arguments{:}) ;
                        self.job_ids(job_index) = this_job_id ;
                        self.has_job_been_submitted(job_index) = true ;
                        is_job_in_progress(job_index) = self.do_actually_submit ;
                    end
                end                                
                has_job_exited = self.has_job_been_submitted & ~is_job_in_progress ;
                exited_job_count = sum(has_job_exited) ;
                newly_exited_job_count = exited_job_count - last_exited_job_count ;
                if do_show_progress_bar ,
                    progress_bar.update(newly_exited_job_count) ;
                end
                have_all_exited = (exited_job_count==job_count) ;
                if ~have_all_exited ,  
                    if self.do_actually_submit ,
                        pause(20) ;
                    end
                    is_time_up = (toc(ticId) > maximum_wait_time) ;
                end
                last_exited_job_count = exited_job_count ;
            end
            job_statuses = get_bsub_job_status(self.job_ids, self.submit_host_name) ;
        end
    end
end
