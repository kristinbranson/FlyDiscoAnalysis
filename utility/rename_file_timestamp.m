function newloc = rename_file_timestamp(filename,sshhost)

if nargin < 2,
  sshhost = '';
end

res = dir(filename);
timestamp = datestr(res.datenum,'yyyymmddTHHMMSS');
newloc = [filename,'.',timestamp];
[success,msg] = mymovefile(filename,newloc,sshhost);
if ~success,
  error('Could not rename existing output file %s: %s',filename,msg);
end

function [success,msg] = mymovefile(filename,newloc,sshhost)

if nargin < 3 || isempty(sshhost),
  [success,msg] = movefile(filename,newloc);
  return;
end

cmd = sprintf('ssh %s ''mv "%s" "%s"''',sshhost,filename,newloc);
[status,result] = unix(cmd);
success = status == 0;
msg = strtrim(result);