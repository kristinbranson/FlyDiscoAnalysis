function saveRegistrationDataToTextFile(file_name, registration_data)
% Saves particular fields from registration_data to the text file named
% file_name.  

fid = fopen(file_name,'w');
cleaner = onCleanup(@()(fclose(fid))) ;

field_names_we_want_to_save = ...
  {'offX','offY','offTheta','scale','bowlMarkerTheta','featureStrengths',...
   'circleCenterX','circleCenterY','circleRadius',...
   'seconds_crop_start','seconds_crop_end','start_frame','end_frame',...
   'flytracker_nnanframes','flytracker_nids0','flytracker_nidsnew'} ;

field_names_we_will_save = intersect(field_names_we_want_to_save, fieldnames(registration_data)) ;
for i = 1:numel(field_names_we_will_save),
  field_name = field_names_we_will_save{i};
  fprintf(fid,'%s,%f\n',field_name,registration_data.(field_name));
end

if isfield(registration_data,'ledIndicatorPoints')
  fprintf(fid,'%s,%d\n','ledX',registration_data.ledIndicatorPoints(1));
  fprintf(fid,'%s,%d\n','ledY',registration_data.ledIndicatorPoints(2));
end
  
end
