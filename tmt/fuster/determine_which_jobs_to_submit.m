function will_submit_from_submittable_index = determine_which_jobs_to_submit(slot_count_from_submittable_index, maximum_slot_count)
    % Determine which of the submittable jobs will be submitted, given how
    % many slots each submittable job needs, and the maximum number of slots we can use.
    submittable_count = length(slot_count_from_submittable_index) ;
    will_submit_from_submittable_index = false(size(slot_count_from_submittable_index)) ;
    slots_used_so_far = 0 ;
    for submittable_index = 1 : submittable_count ,
        slot_count_this_submittable = slot_count_from_submittable_index(submittable_index) ;
        putative_slots_used = slots_used_so_far + slot_count_this_submittable ;
        if putative_slots_used <= maximum_slot_count ,
            will_submit_from_submittable_index(submittable_index) = true ;
            slots_used_so_far = putative_slots_used ;
            if slots_used_so_far >= maximum_slot_count ,
                break
            end
        end                                
    end    
end
