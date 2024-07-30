function ctrax_results_file_from_flytracker_jaaba_trx_file(ctrax_results_mat_name, flytracker_jaaba_trx_file_name)
    flytrack_jaaba_trx_file_as_struct = load(flytracker_jaaba_trx_file_name, '-mat') ;
    ctrax_results_as_struct = ctrax_results_struct_from_flytracker_jaaba_trx_struct(flytrack_jaaba_trx_file_as_struct) ;
    save(ctrax_results_mat_name, '-struct', 'ctrax_results_as_struct') ;
end
