function [success] = ResaveMetadata(metadata,savefilename)
% tested with structcompare

success = true;

% back up metadata file if it exists
% keeps the orginal file's timestamps
if exist(savefilename,'file'),
  bakfilename = [savefilename,'_',datestr(now, 'yyyymmddTHHMMSS'),'.bak'];
  [success1] = copyfile(savefilename,bakfilename,'f');
  if ~success1,
    fprintf('Error backing up metadatafile, aborting ResaveMetadata\n');
    success = false;
    return;
  end
end

if ~exist(bakfilename,'file')
    sprintf('Error backing up metadatafile, aborting ResaveMetadata\n')
    success = false;
    return;
end
  
% open metadata file for writing 
% fails silentely if savefilename is softlink without permissions !

fid = fopen(savefilename,'w');
if fid < 0,
  sprintf('Could not write to experiment metadata file %s',savefilename);
  success = false;
  return;
end

% write the main metadata file
fprintf(fid,'<?xml version="1.0"?>\n');
% name of assay
fprintf(fid,'<experiment assay="%s" ',metadata.assay);
% start datetime
fprintf(fid,'exp_datetime="%s" ',metadata.exp_datetime);
% name of experimenter
fprintf(fid,'experimenter="%s" ',metadata.experimenter);
% lab
fprintf(fid,'lab="%s" ',metadata.lab);
% always same experiment protocol
fprintf(fid,'protocol="%s" ',metadata.protocol);
% screen type
fprintf(fid,'screen_type="%s" ',metadata.screen_type);
% screen reason
fprintf(fid,'screen_reason="%s" ',metadata.screen_reason);
% data capture code version
fprintf(fid,'data_capture_version="%s" ',metadata.data_capture_version);
% optogenetic activation LED protocol name
fprintf(fid,'led_protocol="%s" ',metadata.led_protocol);

fprintf(fid,'>\n');

% session container
fprintf(fid,'  <session id="1">\n');

% apparatus full id and parts
fprintf(fid, '    <apparatus apparatus_id="%s" room="%s" rig="%s" plate="%s" top_plate="%s" visual_surround="%s" bowl="%s" camera="%s" computer="%s" harddrive="%s"/>\n',...
  metadata.apparatus_id,...
  metadata.room,num2str(metadata.rig),num2str(metadata.plate),num2str(metadata.top_plate),num2str(metadata.visual_surround),metadata.bowl,...
  metadata.camera,...
  metadata.computer,...
  metadata.harddrive);
% line name
fprintf(fid,'    <flies line="%s" ', metadata.line);
% effector
fprintf(fid,'effector="%s" ',metadata.effector);
% gender
fprintf(fid,'gender="%s" ',metadata.gender); 
% cross date
fprintf(fid,'cross_date="%s" ',metadata.cross_date);
% flip_date
fprintf(fid,'flip_date="%s" ',metadata.flip_date);

% hours starved
fprintf(fid,'hours_starved="%f" ',metadata.hours_starved);

% barcode
fprintf(fid,'cross_barcode="%d" ',metadata.cross_barcode);
% flip
fprintf(fid,'flip_used="%d" ',metadata.flip_used);
% wish list
fprintf(fid,'wish_list="%d" ',metadata.wish_list);
% robot stock copy. set this to unknown for now
fprintf(fid,'robot_stock_copy="%d" ',metadata.robot_stock_copy);
% count is set to 0 -- won't know this til after tracking
fprintf(fid,'num_flies="0">\n');


% choose rearing protocol based on incubator ID
fprintf(fid,'      <rearing rearing_protocol="%s" ',metadata.rearing_protocol);
fprintf(fid,'rearing_incubator="%s" ',num2str(metadata.rearing_incubator));
fprintf(fid,'/>\n');

% always same handling protocol
fprintf(fid,'      <handling handling_protocol="%s" ',metadata.handling_protocol);
% person who crossed flies
fprintf(fid,'handler_cross="%s" ',metadata.handler_cross);
% person who sorted flies
fprintf(fid,'handler_sorting="%s" ',metadata.handler_sorting);
% time since sorting, in hours
fprintf(fid,'hours_sorted="%f" ',metadata.hours_sorted);

% person who moved flies to starvation material
fprintf(fid,'handler_starvation="%s" ',metadata.handler_starvation);

% seconds between bringing vials into hot temperature environment and
% experiment start
fprintf(fid,'seconds_shiftflytemp="%f" ',metadata.seconds_shiftflytemp);
% seconds between loading flies into arena and experiment start
fprintf(fid,'seconds_fliesloaded="%f" ',metadata.seconds_fliesloaded);
% number of observed dead flies
fprintf(fid,'num_flies_dead="%d" ',metadata.num_flies_dead);
% number of observed damaged flies
fprintf(fid,'num_flies_damaged="%d" ',metadata.num_flies_damaged);
fprintf(fid,'/>\n');
fprintf(fid,'    </flies>\n');
fprintf(fid,'  </session>\n');
% temperature and humidity measured from precon sensor
fprintf(fid,'  <environment temperature="%f" ',metadata.temperature);
fprintf(fid,'humidity="%f" />\n',metadata.humidity);
% notes entered
% deal with multi-line notes
fprintf(fid,'  <notes_behavioral>%s</notes_behavioral>\n',metadata.notes_behavioral);
% % deal with multi-line notes
% if iscell(handles.TechnicalNotes),
%   TechnicalNotes = handles.TechnicalNotes;
% else
%   TechnicalNotes = cellstr(handles.TechnicalNotes);
% end
% TechnicalNotes = sprintf('%s\\n',TechnicalNotes{:});
% TechnicalNotes = TechnicalNotes(1:end-2);
fprintf(fid,'  <notes_technical>%s</notes_technical>\n',metadata.notes_technical);
fprintf(fid,'  <notes_keyword></notes_keyword>\n');
% flags entered
fprintf(fid,'  <flag_review>%d</flag_review>\n',metadata.flag_review);
fprintf(fid,'  <flag_redo>%d</flag_redo>\n',metadata.flag_redo);
fprintf(fid,'  <flag_aborted>%d</flag_aborted>\n',metadata.flag_aborted);
fprintf(fid,'  <flag_legacy>0</flag_legacy>\n');

fprintf(fid,'</experiment>\n');

fclose(fid);






