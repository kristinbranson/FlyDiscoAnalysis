%%

rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
filestr = 'movie.ufmf';
datetime_format = 'yyyymmddTHHMMSS';
max_dt = 7;
min_dt = 0;
outfilename = 'TimestampVsExpDateTimeCheckResults20120105.txt';


%%

files = dir(rootdatadir);

%%

fid = fopen(outfilename,'w');

for i = 1:numel(files),
  exp_datetime = regexp(files(i).name,'_(\d{8}T\d{6})$','once','tokens');
  if isempty(exp_datetime),
    fprintf('Skipping %s, not an experiment directory.\n',files(i).name);
    continue;
  end
  exp_datetime = exp_datetime{1};
  fn = fullfile(rootdatadir,files(i).name,filestr);
  if ~exist(fn,'file'),
    fprintf('Skipping %s, does not contain file %s\n',files(i).name,filestr);
    continue;
  end
  exp_datenum = datenum(exp_datetime,datetime_format);
  if isnan(exp_datenum),
    warning('Could not parse exp_datetime %s, skipping\n',exp_datetime);
    continue;
  end
  moviefile = dir(fn);
  dt = moviefile.datenum - exp_datenum;
  if dt > max_dt || dt < min_dt,
    
    % check if this is on archive
    [~,link] = unix(sprintf('readlink %s',fullfile(rootdatadir,files(i).name)));
    if ~isempty(link) && ~isempty(strfind(link,'archive')),
      fprintf('Skipping %s, archived\n',files(i).name);
      continue;
    end
    
    fprintf(fid,'%s: dt = %.1f days, exp_datetime = %s, file timestamp = %s, ',...
      files(i).name,dt,exp_datetime,datestr(moviefile.datenum,datetime_format));
    same_datetime_files = dir(fullfile(rootdatadir,sprintf('*_%s',exp_datetime)));
    same_datetime_fns = setdiff({same_datetime_files.name},{files(i).name});
    if isempty(same_datetime_fns),
      fprintf(fid,'No directories with same exp_datetime\n');
    else
      s = sprintf('%s ',same_datetime_fns{:});
      fprintf(fid,'possible matches: %s\n',s);
    end
    
  end
end

fclose(fid);

%% find experiments with the same datetime as problem experiments

suspexps = {};
fid = fopen(outfilename,'r');
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  i = find(s==':',1);
  if isempty(i),
    continue;
  end
  expname = s(1:i-1);
  
  exp_datetime = regexp(files(i).name,'_(\d{8}T\d{6})$','once','tokens');
  exp_datetime = exp_datetime{1};
  
  dir(fullfile(rootdatadir,sprintf('*_%s',exp_datetime)));
  
end
fclose(fid);