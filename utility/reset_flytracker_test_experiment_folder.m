function reset_flytracker_test_experiment_folder()
  template_test_experiment_folder_path = '/groups/branson/bransonlab/adam/experiments/emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2_short_with_ft_template' ;
  test_experiment_folder_path = '/groups/branson/bransonlab/adam/experiments/emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2_short_with_ft' ;
  if exist(test_experiment_folder_path, 'file') ,
    rmdir(test_experiment_folder_path, 's') ;
  end
  mkdir(test_experiment_folder_path) ;
  copyfile(template_test_experiment_folder_path, test_experiment_folder_path) ;  
end
