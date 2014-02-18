function [filenamecurr] = ShowDNLineInfo(line_names,isdnexpr,dnnames,varargin)

persistent metadata;
timestamp = datestr(now,'yyyymmddTHHMMSS');

[metadata,lineresultsdir,outlineresultsdir,groupname,filenamecurr,extraimnames,extraimheights,olympiaddir,outolympiaddir] = myparse(varargin,...
  'metadata',metadata,...
  'lineresultsdir','/groups/branson/bransonlab/projects/olympiad/LineResults',...
  'outlineresultsdir','/groups/branson/bransonlab/projects/olympiad/LineResults',...
  'groupname',sprintf('Group selected %s',timestamp),...
  'outfilename','',...
  'extraimnames',{},...
  'extraimheights',[],...
  'olympiaddir','',...
  'outolympiaddir','');

if isempty(metadata),
  [f,p] = uigetfile('*.mat','Select metadata mat file');
  if ~ischar(f),
    return;
  end
  metadatafile = fullfile(p,f);
  metadata = load(metadatafile,'metadata');
  metadata = metadata.metadata;
end

nlines = numel(line_names);

% create an html page with links to all experiments
if isempty(filenamecurr),
  filenamecurr = fullfile(tempdir,sprintf('selectedgroup_%s.html',timestamp));
end

fid = fopen(filenamecurr,'w');
if fid < 1,
  error('Could not open temporary file %s for writing.\n');
end

fprintf(fid,'<html>\n<title>%s</title>\n<body>\n',groupname);
fprintf(fid,'<head>\n');
fprintf(fid,'<style>\n');
fprintf(fid,'table\n');
fprintf(fid,'{\n');
fprintf(fid,'border-collapse:collapse;\n');
fprintf(fid,'}\n');
fprintf(fid,'table, td, th\n');
fprintf(fid,'{\n');
fprintf(fid,'border:1px solid black;\n');
fprintf(fid,'}\n');
fprintf(fid,'</style>\n');
fprintf(fid,'</head>\n');

for i = 1:numel(extraimnames),
  if i <= numel(extraimheights),
    fprintf(fid,'<p><a href="%s"><img src="%s" height="%d"></a></p>\n',...
      extraimnames{i},extraimnames{i},extraimheights(i));
  else
    fprintf(fid,'<p><a href="%s"><img src="%s"></a></p>\n',...
      extraimnames{i},extraimnames{i});
  end
end
        
for i = 1:nlines,
  fprintf(fid,'<h1>%s</h1>\n',line_names{i});
  
  % show line results
  resultsdirname = fullfile(lineresultsdir,line_names{i});
  outresultsdirname = fullfile(outlineresultsdir,line_names{i});
  if ~exist(resultsdirname,'dir'),
    fprintf(fid,'<p>Line results not found.</p>\n');
  else
    imname = fullfile(outresultsdirname,'stats_basic.png');
    fprintf(fid,'<p><a href="%s"><img src="%s" height="300"></a></p>\n',...
      imname,imname);
    fprintf(fid,'<p><a href="%s">Line behavior plots</a></p>\n',outresultsdirname);
  end
  
  idx = find(strcmp({metadata.line_name},line_names{i}));
  if isempty(idx),
    fprintf(fid,'<p>No experiments found.</p>\n');
  else
    fprintf(fid,'<ul>\n');
    for j = idx(:)',
      [~,name] = fileparts(metadata(j).file_system_path);
      if isempty(olympiaddir),
        expdir = metadata(j).file_system_path;
      else
        expdir = fullfile(olympiaddir,name);
      end
      if isempty(outolympiaddir),
        outexpdir = metadata(j).file_system_path;
      else
        outexpdir = fullfile(outolympiaddir,'fly_bowl','bowl_data',name);
      end
      
      moviename = fullfile(expdir,...
        sprintf('ctrax_results_movie_%s.avi',name));
      outmoviename = fullfile(outexpdir,...
        sprintf('ctrax_results_movie_%s.avi',name));
        
      plotsname = fullfile(expdir,'analysis_plots');
      outplotsname = fullfile(outexpdir,'analysis_plots');
      if ~exist(moviename,'file') && ~exist(plotsname,'file'),
        continue;
      end
      fprintf(fid,'  <li>%s: ',name);
      if exist(moviename,'file'),
        fprintf(fid,'<a href="file://%s">Ctrax results movie</a>',outmoviename);
      end
      if exist(moviename,'file') && exist(plotsname,'dir'),
        fprintf(fid,', ');
      end
      if exist(plotsname,'dir'),
        fprintf(fid,'<a href="file://%s">Analysis plots</a>',outplotsname);
      end
      fprintf(fid,'</li>\n');
      
    end
    fprintf(fid,'</ul>\n');
  end
end
      
fprintf(fid,'<table>\n');
fprintf(fid,'<tr><th>line</th>');
for i = 1:numel(dnnames),
  fprintf(fid,'<th>%s</th> ',dnnames{i});
end
fprintf(fid,'</tr>\n');
for i = 1:nlines,
  fprintf(fid,'<td>%s</td>',line_names{i});
  for j = 1:numel(dnnames),
    fprintf(fid,'<td>%d</td>',isdnexpr(i,j));
  end
  fprintf(fid,'</tr>\n');
end
fprintf(fid,'</table>\n');
    
fprintf(fid,'</body>\n</html>\n');
fclose(fid);

if ~exist(filenamecurr,'file'),
  warning('Could not open temporary file %s',filenamecurr);
else
  % open this page
  web(filenamecurr,'-browser');
end
