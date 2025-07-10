function MakeLineImageDirectory(line_name,linestats,setstats,expstats,varargin)

[rootdir,deleteold,lineresultsrootdir,...
  statscsvfilestr,plotdirstr,...
  statfns,distweights,...
  anatomyurlprefix,statfns_table] = myparse(varargin,...
  'rootdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlResults',...
  'deleteold',true,...
  'lineresultsrootdir','/groups/branson/bransonlab/projects/olympiad/LineResults',...
  'statscsvfilestr','data.csv',...
  'plotdirstr','analysis_plots',...
  'statfns',{},...
  'distweights',[],...
  'anatomyurlprefix','http://flweb.janelia.org/cgi-bin/view_flew_imagery.cgi?line=',...
  'statfns_table',MakeLineImage_SetStats());

if isempty(statfns),
  statfns = fieldnames(expstats.zscores);
end
nstats = numel(statfns);

linedir = fullfile(rootdir,line_name);
if exist(linedir,'dir') && deleteold,
  fprintf('Deleting exist line directory %s...\n',linedir);
  unix(sprintf('rm -r %s',linedir));
end
if ~exist(linedir,'dir'),
  mkdir(linedir);
end

linei = find(strcmp({linestats.metadata.line_name},line_name));
if isempty(linei),
  error('Could not find line data for line %s',line_name);
end

% create a csv file with the per-frame stats
csvfilename = fullfile(linedir,statscsvfilestr);
[success,msg] = SaveStatsCSVFile(linestats,linei,csvfilename,'statfns',statfns);
if ~success,
  error(msg);
end

% copy over the plots
lineresultsdir = fullfile(lineresultsrootdir,line_name);
if ~exist(lineresultsdir,'dir'),
  error('Line results directory %s does not exist',lineresultsrootdir);
end

% create the output plot directory
lineplotdir = fullfile(linedir,plotdirstr);
mkdir(lineplotdir);

unix(sprintf('cp %s/*.png %s/.',lineresultsdir,lineplotdir));

% create the experiment directories
expis = find(strcmp({expstats.metadata.line_name},line_name));
nexps = numel(expis);

for ii = 1:nexps,
  
  i = expis(ii);
  [~,experiment_name] = fileparts(expstats.metadata(i).file_system_path);
  expdir = fullfile(linedir,experiment_name);
  
  % create the directory
  mkdir(expdir);
  
  % save the stats file
  csvfilename = fullfile(expdir,statscsvfilestr);
  [success,msg] = SaveStatsCSVFile(expstats,i,csvfilename,'statfns',statfns);
  if ~success,
    error(msg);
  end
  
  % copy over the plots
  inexpdir = fullfile(expstats.metadata(i).file_system_path);
  if ~exist(inexpdir,'dir'),
    error('Experiment directory %s does not exist',inexpdir);
  end
  inplotdir = fullfile(inexpdir,'analysis_plots');
  if ~exist(inplotdir,'dir'),
    error('analysis_plots directory does not exist for experiment directory %s',inexpdir);
  end
  plotdir = fullfile(expdir,plotdirstr);
  mkdir(plotdir);
  unix(sprintf('cp %s/*.png %s/.',inplotdir,plotdir));
  
  % copy the metadata file
  inmetadatafile = fullfile(inexpdir,'Metadata.xml');
  if ~exist(inmetadatafile,'file'),
    error('Metadata file %s does not exist',inmetadatafile);
  end
  unix(sprintf('cp %s %s/.',inmetadatafile,expdir));
  
  % copy the results video
  inresultsvideo = fullfile(inexpdir,sprintf('ctrax_results_movie_%s.avi',experiment_name));
  if ~exist(inresultsvideo,'file'),
    warning('Results video %s does not exist',inresultsvideo);
  else
    unix(sprintf('cp %s %s/.',inresultsvideo,expdir));
  end

end

%% create the html page
htmlfile = fullfile(linedir,'index.html');
fid = fopen(htmlfile,'w');
if fid < 1,
  error('Could not open file %s for writing',htmlfile);
end

% output the header

fprintf(fid,'<html>\n');
fprintf(fid,'<head>\n');
fprintf(fid,'<title>%s</title>\n',line_name);
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n');
fprintf(fid,'<link rel="stylesheet" type="text/css" charset="utf-8" media="all" href="../styles/common.css">\n');
fprintf(fid,'<link rel="stylesheet" type="text/css" charset="utf-8" media="screen" href="../styles/screen.css">\n');
fprintf(fid,'<link rel="stylesheet" type="text/css" charset="utf-8" media="print" href="../styles/print.css">\n');
fprintf(fid,'<script src="../js/jquery-1.10.2.min.js"></script>\n');
fprintf(fid,'<script src="../js/jquery.tablesorter.js"></script>\n');
fprintf(fid,'<script src="../js/jquery.tablesorter.pager.js"></script>\n');

fprintf(fid,'<script type="text/javascript" id="js">$(document).ready(function() {\n');
fprintf(fid,'// call the tablesorter plugin\n');
fprintf(fid,'$("table").tablesorter({\n');
fprintf(fid,'// sort on the third column, order desc\n');
fprintf(fid,'sortList: [[2,1]]\n');
fprintf(fid,'});\n');
fprintf(fid,'}); </script>\n');


fprintf(fid,'</head>\n');
fprintf(fid,'<body>\n');

fprintf(fid,'<h2>%s: Fly Bowl neural activation screen behavior phenotypes</h2>\n',line_name);

% basic stats plot
fprintf(fid,'<div id="stats_basic_plot">\n');
fprintf(fid,'<p align="center"><a href="%s"><img src="%s" width="500px" style="float:center"/></a></p>\n',...
  fullfile(plotdirstr,'stats_basic.png'),...
  fullfile(plotdirstr,'stats_basic.png'));
fprintf(fid,'</div>\n');

fprintf(fid,'<div id="links">\n');
fprintf(fid,'<h3>Links</h3>\n');
fprintf(fid,'<p><ul>\n');

% link to results video

% choose the most "normal" experiment to use as the sample results video

if isempty(distweights),
  distweights = ones(1,nstats);
end

expdata = nan(nexps,nstats);
for i = 1:nstats,
  fn = statfns{i};
  expdata(:,i) = expstats.zscores.(fn)(expis);
end

linedata = nan(1,nstats);
for i = 1:nstats,
  fn = statfns{i};
  linedata(i) = linestats.zscores.(fn)(linei);
end

d = nanmean(bsxfun(@times,abs(bsxfun(@minus,expdata,linedata)),distweights),2);
[~,i] = min(d);
expiclosest = expis(i);

expday = expstats.metadata(expiclosest).exp_datetime(1:8);
fprintf(fid,'<li><a href="%s">Results video</a> for experiment %s%s</li>\n',...
  fullfile(experiment_name,sprintf('ctrax_results_movie_%s.avi',experiment_name)),...
  expday,expstats.metadata(expiclosest).bowl);
fprintf(fid,'<li><a href="%s">Analysis plots</a></li>\n',plotdirstr);

shortlinename = regexprep(line_name,'^GMR_','R');
shortlinename = regexprep(shortlinename,'_A._01$','');

fprintf(fid,'<li><a href="%s%s">Fly Light imagery</a></li>\n',anatomyurlprefix,shortlinename);
fprintf(fid,'<li><a href="%s">Data csv file</a></li>\n',statscsvfilestr);

fprintf(fid,'</ul></p>\n');
fprintf(fid,'</div>\n');

fprintf(fid,'<div id="experiments">\n');

fprintf(fid,'<p><h3>Individual trials</h3></p>\n');
fprintf(fid,'<p><ul>\n');

setis = find(strcmp({setstats.metadata.line_name},line_name));
nsets = numel(setis);

for ii = 1:nsets,
  seti = setis(ii);
  setday = setstats.metadata(seti).exp_datetime(1:8);  
  if nsets > 1,
    fprintf(fid,'%s\n<ul>\n',setday);
  end
  
  expis1 = find(strcmp({expstats.metadata.set},setstats.metadata(seti).set));
  [~,order] = sort({expstats.metadata(expis1).bowl});
  expis1 = expis1(order);
  for jj = 1:numel(expis1),
    expj = expis1(jj);
    [~,experiment_name] = fileparts(expstats.metadata(expj).file_system_path);
    fprintf(fid,'<li>%s: <a href="%s">Results video</a>&nbsp&nbsp|&nbsp&nbsp',...
      expstats.metadata(expj).bowl,...
      fullfile(experiment_name,sprintf('ctrax_results_movie_%s.avi',experiment_name)));
    fprintf(fid,'<a href="%s">Summary statistics</a>&nbsp&nbsp|&nbsp&nbsp',...
      fullfile(experiment_name,plotdirstr,'stats_basic.png'));
    fprintf(fid,'<a href="%s">Other plots</a>&nbsp&nbsp|&nbsp&nbsp',fullfile(experiment_name,plotdirstr));
    fprintf(fid,'<a href="%s">Data csv file</a>&nbsp&nbsp|&nbsp&nbsp',fullfile(experiment_name,statscsvfilestr));
    fprintf(fid,'<a href="%s">Metadata file</a>',fullfile(experiment_name,'Metadata.xml'));
    fprintf(fid,'</li>\n');
  end
  
  if nsets > 1,
    fprintf(fid,'</ul>\n');
  end
end  


fprintf(fid,'</ul></p>\n');

fprintf(fid,'</div>\n');


% table of p-values
fprintf(fid,'<div id="resultstable">\n');

fprintf(fid,'<h3>Behavior statistics</h3>\n');

fprintf(fid,'<table cellspacing="1" class="tablesorter">\n');
fprintf(fid,'<thead>\n');
fprintf(fid,'<tr><th>Behavior statistic</th><th>Raw value</th><th>Z-score</th><th>p-value: &gt; control</th><th>p-value: &lt; control</th></tr>\n');
fprintf(fid,'</thead>\n');
fprintf(fid,'<tbody>\n');

statfns_table = intersect(statfns_table,statfns);
for stati = 1:numel(statfns_table),
  fn = statfns_table{stati};
  if isnan(linestats.means.(fn)(linei)),
    continue;
  end
  fprintf(fid,'<tr><td>%s</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td></tr>\n',fn,linestats.means.(fn)(linei),...
    linestats.zscores.(fn)(linei),min(1,linestats.pvalue_bigger.(fn)(linei)),...
    min(1,linestats.pvalue_smaller.(fn)(linei)));
end
fprintf(fid,'</tbody>\n');
fprintf(fid,'</table>\n');

fprintf(fid,'</div>\n');



fprintf(fid,'</body>\n');
fprintf(fid,'</html>\n');

fclose(fid);