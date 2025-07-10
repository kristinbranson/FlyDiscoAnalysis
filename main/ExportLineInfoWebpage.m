function [hfigs,hax,filenamecurr] = ExportLineInfoWebpage(line_names,outdir,varargin)

persistent metadata;
anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
lineresultsdir = '/groups/branson/bransonlab/projects/olympiad/LineResults';
timestamp = datestr(now,'yyyymmddTHHMMSS');
hfigs_default = [];
persistent imdata;

[metadata,imdata,showanatomy,showbehavior,showaverageanatomy,hfigs1,anatomydir,lineresultsdir,...
  groupname,int_manual,filenamecurr] = myparse(varargin,...
  'metadata',metadata,...
  'imdata',imdata,...
  'showanatomy',true,...
  'showbehavior',true,...
  'showaverageanatomy',false,...
  'hfigs',hfigs_default,...
  'anatomydir',anatomydir,...
  'lineresultsdir',lineresultsdir,...
  'groupname',sprintf('Group selected %s',timestamp),...
  'int_manual',[],...
  'outfilename','');
if numel(hfigs1) >= numel(hfigs_default),
  hfigs = hfigs1;
elseif isempty(hfigs1),
  hfigs = hfigs_default;
else
  hfigs = hfigs_default;
  hfigs(1:numel(hfigs1)) = hfigs1;
end

copyflags = 'nL';

% create output directories
if ~exist(outdir,'dir'),
  mkdir(outdir);
end
nlines = numel(line_names);
for i = 1:nlines,
  outlinedircurr = fullfile(outdir,line_names{i});
  if ~exist(outlinedircurr,'dir'),
    mkdir(outlinedircurr);
  end
end

if isempty(metadata),
  [f,p] = uigetfile('*.mat','Select metadata mat file');
  if ~ischar(f),
    return;
  end
  metadatafile = fullfile(p,f);
  metadata = load(metadatafile,'metadata');
  metadata = metadata.metadata;
end
if isempty(imdata),
  [f,p] = uigetfile('*.mat','Select imdata mat file','ImageryData20130824.mat');
  if ~ischar(f),
    return;
  end
  imagerydatafile = fullfile(p,f);
  imdata = load(imagerydatafile,'imdata');
  imdata = imdata.imdata;
end  

if isempty(filenamecurr),
  filenamecurr = sprintf('selectedgroup_%s.html',timestamp);
end
fid = fopen(fullfile(outdir,filenamecurr),'w');
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

% create an html page with links to all experiments
if showbehavior || (showanatomy && ~isempty(int_manual)),
  
  
  if showbehavior,
    for i = 1:nlines,
      fprintf(fid,'<h1>%s</h1>\n',line_names{i});
      lineresultsdircurr = fullfile(lineresultsdir,line_names{i});
      if exist(lineresultsdircurr,'dir'),
        outlineresultsdircurr = fullfile(line_names{i},'LineResultPlots');
        if ~exist(fullfile(outdir,outlineresultsdircurr),'file'),
          mkdir(fullfile(outdir,outlineresultsdircurr));
        end
        cmd = sprintf('cp -%s %s/stat*png %s/.',copyflags,lineresultsdircurr,fullfile(outdir,outlineresultsdircurr));
        fprintf('%s\n',cmd);
        unix(cmd);
        lineresultsimfile = fullfile(outlineresultsdircurr,'stats_basic.png');
        fprintf(fid,'<p><a href="%s"><img src="%s" height="500"></a></p>\n',lineresultsimfile,lineresultsimfile);
        fprintf(fid,'<p><a href="%s">Line result plots</a></p>\n',outlineresultsdircurr);
      else
        fprintf(fid,'<p>Line results plots not found.</p>\n');
      end
      
      idx = find(strcmp({metadata.line_name},line_names{i}));
      if isempty(idx),
        fprintf(fid,'<p>No experiments found.</p>\n');
      else
        fprintf(fid,'<ul>\n');
        for j = idx(:)',
          [~,name] = fileparts(metadata(j).file_system_path);
          outexpdircurr = fullfile(line_names{i},name);
          if ~exist(fullfile(outdir,outexpdircurr),'dir')
            mkdir(fullfile(outdir,outexpdircurr));
          end
          moviename = fullfile(metadata(j).file_system_path,...
            sprintf('ctrax_results_movie_%s.avi',name));
          outmoviename = fullfile(outexpdircurr,sprintf('ctrax_results_movie_%s.avi',name));
          if exist(moviename,'file'),
            cmd = sprintf('cp -%s %s %s/.',copyflags,moviename,fullfile(outdir,outexpdircurr));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          plotsname = fullfile(metadata(j).file_system_path,'analysis_plots');
          outplotsname = fullfile(outexpdircurr,'analysis_plots');
          if exist(plotsname,'dir'),
            if ~exist(fullfile(outdir,outplotsname),'dir'),
              mkdir(fullfile(outdir,outplotsname));
            end
            cmd = sprintf('cp -%s %s/stat*.png %s/.',copyflags,plotsname,fullfile(outdir,outplotsname));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          if ~exist(moviename,'file') && ~exist(plotsname,'file'),
            continue;
          end
          fprintf(fid,'  <li>%s: ',name);
          if exist(moviename,'file'),
            fprintf(fid,'<a href="%s">Ctrax results movie</a>',outmoviename);
          end
          if exist(moviename,'file') && exist(plotsname,'dir'),
            fprintf(fid,', ');
          end
          if exist(plotsname,'dir'),
            fprintf(fid,'<a href="%s">Analysis plots</a>',outplotsname);
          end
          fprintf(fid,'</li>\n');
          
        end
        fprintf(fid,'</ul>\n');
      end
      idxanat = find(strcmp({imdata.line},line_names{i}));
      if isempty(idxanat),
        fprintf(fid,'<p>No images found.</p>\n');
      else
        fprintf(fid,'<ul>\n');
        outanatdir = fullfile(line_names{i},'FlyLightImages');
        if ~exist(fullfile(outdir,outanatdir),'dir'),
          mkdir(fullfile(outdir,outanatdir));
        end
        for j = idxanat(:)',
          name = imdata(j).name;
          maxproj_file = imdata(j).maxproj_file_system_path;
          maxproj_ch2_file = regexprep(maxproj_file,'_total.jpg$','_ch2_total.jpg');
          reg_file = regexprep(maxproj_file,'_total.jpg$','.reg.local.jpg');
          translation_file = imdata(j).translation_file_path;
          if ~exist(maxproj_file,'file') && ~exist(translation_file,'file'),
            continue;
          end
          [~,tmpname] = myfileparts(maxproj_file);
          outmaxproj_file = fullfile(outanatdir,tmpname);
          if exist(maxproj_file,'file'),
            cmd = sprintf('cp -%s %s %s/.',copyflags,maxproj_file,fullfile(outdir,outanatdir));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          [~,tmpname] = myfileparts(maxproj_ch2_file);
          outmaxproj_ch2_file = fullfile(outanatdir,tmpname);
          if exist(maxproj_ch2_file,'file'),
            cmd = sprintf('cp -%s %s %s/.',copyflags,maxproj_ch2_file,fullfile(outdir,outanatdir));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          [~,tmpname] = myfileparts(reg_file);
          outreg_file = fullfile(outanatdir,tmpname);
          if exist(reg_file,'file'),
            cmd = sprintf('cp -%s %s %s/.',copyflags,reg_file,fullfile(outdir,outanatdir));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          [~,tmpname] = myfileparts(translation_file);
          outtranslation_file = fullfile(outanatdir,tmpname);
          if exist(translation_file,'file'),
            cmd = sprintf('cp -%s %s %s/.',copyflags,translation_file,fullfile(outdir,outanatdir));
            fprintf('%s\n',cmd);
            unix(cmd);
          end
          
          fprintf(fid,'  <li>%s: ',name);
          if exist(maxproj_file,'file'),
            fprintf(fid,'<a href="%s">Max projection image</a>',outmaxproj_file);
          end
          if exist(maxproj_ch2_file,'file') && ~strcmp(maxproj_file,maxproj_ch2_file),
            fprintf(fid,', <a href="%s">Max proj, channel 2</a>',outmaxproj_ch2_file);
          end
          if exist(reg_file,'file') && ~strcmp(maxproj_file,reg_file),
            fprintf(fid,', <a href="%s">Registered image</a>',outreg_file);
          end
          if exist(translation_file,'file'),
            fprintf(fid,', <a href="%s">Translation video</a>',outtranslation_file);
          end
          fprintf(fid,'</li>\n');
          
        end
        fprintf(fid,'</ul>\n');
      end
    end
  end
  if showanatomy && ~isempty(int_manual),
    
    compartments = setdiff(fieldnames(int_manual),{'line_names'});
    
    fprintf(fid,'<table>\n');
    fprintf(fid,'<tr><th>line</th>');
    for i = 1:numel(compartments),
      fprintf(fid,'<th>%s</th> ',compartments{i});
    end
    fprintf(fid,'</tr>\n');
    [ism,idx] = ismember(line_names,int_manual.line_names);
    for i = 1:nlines,
      fprintf(fid,'<td>%s</td>',line_names{i});
      if ism(i),
        for j = 1:numel(compartments),
          fprintf(fid,'<td>%.1f</td>',int_manual.(compartments{j})(idx(i)));
        end
      else
        for j = 1:numel(compartments),
          fprintf(fid,'<td>%d</td>',nan);
        end
      end
      fprintf(fid,'</tr>\n');
    end
    fprintf(fid,'</table>\n');
  end
  
end

hax = [];
if showanatomy,
  
  if numel(hfigs) < 1,
    hfigs(1) = figure;
  end
  figure(hfigs(1));
  clf;
  set(hfigs(1),'Units','pixels','Position',[10 10 1800 1200]);
  nc = ceil(sqrt(nlines));
  nr = ceil((nlines)/nc);
  hax = createsubplots(nr,nc,[.01,.01],hfigs(1));
  if numel(hax) > nlines,
    delete(hax(nlines+1:end));
    hax = hax(1:nlines);
  end

  hwait = mywaitbar(0,'Reading max projection images');
  nlinesplotted = 0;
  didplot = false(1,nlines);
  for i = 1:nlines,
    filename = fullfile(anatomydir,sprintf('meanim_%s.png',line_names{i}));
    hwait = mywaitbar((i-1)/(nlines+1),hwait,sprintf('Reading %s...\n',line_names{i}));
    if exist(filename,'file'),
      im = imread(filename);
      image(im,'Parent',hax(i));
      axis(hax(i),'image');
      nlinesplotted = nlinesplotted+1;
      didplot(i) = true;
    else
      fprintf('Line %s has no image\n',line_names{i});
    end
    title(hax(i),line_names{i},'Interpreter','none');
    set(hax(i),'XTick',[],'YTick',[]);
    drawnow;
  end
  
  [~,name] = fileparts(filenamecurr);
  anatfilename = sprintf('%s_anatomyimages.png',name);
  savefig(fullfile(outdir,anatfilename),hfigs(1),'png');
  fprintf(fid,'<h1>Anatomy images</h1>\n');
  fprintf(fid,'<p><a href="%s"><img src="%s" width="1000"></a></p>\n',anatfilename,anatfilename);

  if nlinesplotted > 0,
    impixelinfo(hfigs(1));
    linkaxes(hax(didplot));
  end
  
  if ishandle(hwait),
    delete(hwait);
  end
  
end

if showaverageanatomy,
  mu = single(0);
  nread = 0;
  hwait = mywaitbar(0,'Reading max projection images');
  for i = 1:nlines,
    hwait = mywaitbar((i-1)/(nlines+1),hwait,sprintf('Reading %s...\n',line_names{i}));
    filename = fullfile(anatomydir,sprintf('meanim_%s.mat',line_names{i}));
    if exist(filename,'file'),
      tmp = load(filename,'meanim');
      mu = mu + tmp.meanim;
      nread = nread+1;
    end
  end
  mu = mu / nread;
  if numel(hfigs) < 2,
    hfigs(2) = figure;
  end

  figure(hfigs(2));
  clf;
  set(hfigs(2),'Units','pixels','Position',[10,10,1400 700]);
  hax(end+1) = gca;
  imagesc(max(mu,[],3)','Parent',hax(nlines+1));
  axis(hax(nlines+1),'image');
  title(hax(nlines+1),'Mean image','Interpreter','none');
  colormap(hfigs(2),kjetsmooth(256));
  set(hax(nlines+1),'XTick',[],'YTick',[]);
  colorbar('peer',hax(nlines+1));
  impixelinfo(hfigs(2));
  linkaxes(hax);
  
  [~,name] = fileparts(filenamecurr);
  anatfilename = sprintf('%s_averageanatomyimage.png',name);
  savefig(fullfile(outdir,anatfilename),hfigs(1),'png');
  fprintf(fid,'<h1>Average anatomy image</h1>\n');
  fprintf(fid,'<p><a href="%s"><img src="%s" width="1000"></a></p>\n',anatfilename,anatfilename);

  
  if ishandle(hwait),
    delete(hwait);
  end

end

fprintf(fid,'</body>\n</html>\n');
fclose(fid);
if ~exist(fullfile(outdir,filenamecurr),'file'),
  warning('Could not create temporary file %s',fullfile(outdir,filenamecurr));
else
  % open this page
  web(['file://',filenamecurr],'-browser');
end
