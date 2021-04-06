input_folder_path = ...
    fullfile('/groups/branson/bransonlab/taylora/flydisco/experiments', ...
             'emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2') ;
output_folder_path = 'emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2_shortened' ;

if ~exist(output_folder_path, 'file') ,
  mkdir(output_folder_path) ;
end

file_names_to_copy = {'metaData.xml', 'protocol.mat'} ;
for i = 1 : length(file_names_to_copy) ,
  file_name = file_names_to_copy{i} ;
  source_path = fullfile(input_folder_path, file_name) ;
  target_path = fullfile(output_folder_path, file_name) ;
  if exist(target_path, 'file') ,
    delete(target_path) ;
  end
  copyfile(source_path, target_path) ;
end

video_file_name = 'movie.ufmf' ;
source_video_path = fullfile(input_folder_path, video_file_name) ;
target_video_path = fullfile(output_folder_path, video_file_name) ;

if exist(target_video_path, 'file') ,
  delete(target_video_path) ;
end

desired_frame_count = 2000 ;
shorten_ufmf(source_video_path, target_video_path, desired_frame_count) ;



% try to read frames from the resulting file
header = ufmf_read_header(target_video_path) ;
