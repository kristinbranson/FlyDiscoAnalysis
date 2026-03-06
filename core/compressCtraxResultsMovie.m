function compressCtraxResultsMovie(expdir, avifilestr, avi_file_path, height, width, subtitlefile)

% Tweak the frame size
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);

% Generate the two command-line commands we need
mp4_file_path = fullfile(expdir, [avifilestr,'.mp4']) ;
nowstr = datestr(now(),'yyyymmddTHHMMSSFFF');
passlogfile = sprintf('%s_%s',avi_file_path,nowstr);
ffmpeg_command = 'env -u LD_LIBRARY_PATH /usr/bin/ffmpeg' ;  
cmd1 = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k -vf "subtitles=%s:force_style=''FontSize=10,FontName=Helvetica''" -pass 1 -f mp4 /dev/null',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlefile);
cmd2 = sprintf('%s -i %s -y -passlogfile %s -c:v h264 -pix_fmt yuv420p -s %dx%d -b:v 1600k -vf "subtitles=%s:force_style=''FontSize=10,FontName=Helvetica''" -pass 2 -f mp4 %s',...
  ffmpeg_command, avi_file_path,passlogfile,newwidth,newheight,subtitlefile,mp4_file_path);

% Run the first command, deal with any error
status1 = system(cmd1);
if status1 ~= 0,
  fprintf('*****\n');
  warning('ffmpeg first pass failed.');
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd1);
  fprintf('then\n');
  fprintf('%s\n',cmd2);
  fprintf('then delete %s %s* %s\n',avi_file_path,passlogfile,subtitlefile);
  fprintf('*****\n');
  return
end

% Run the second command, deal with any error
status2 = system(cmd2);
if status2 ~= 0,
  fprintf('*****\n');
  warning('ffmpeg second pass failed.');
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd2);
  fprintf('then delete %s %s* %s\n',avi_file_path,passlogfile,subtitlefile);
  fprintf('*****\n');
  return
end

% If both commands succeeded, clean up the intermediate files that are no longer needed
delete(avi_file_path);
delete(subtitlefile);
fname = [passlogfile '-*.log'];
if isscalar(dir(fname))
  delete(fname);
end
fname = [passlogfile '-*.log.mbtree'];
if isscalar(dir(fname))
  delete(fname);
end
