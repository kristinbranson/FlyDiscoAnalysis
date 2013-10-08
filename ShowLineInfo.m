function [hfigs,hax,filenamecurr] = ShowLineInfo(line_names,varargin)

persistent metadata;
anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
timestamp = datestr(now,'yyyymmddTHHMMSS');
hfigs_default = [];
persistent imdata;

[metadata,imdata,showanatomy,showbehavior,showaverageanatomy,hfigs1,anatomydir,...
  groupname,int_manual] = myparse(varargin,...
  'metadata',metadata,...
  'imdata',imdata,...
  'showanatomy',true,...
  'showbehavior',true,...
  'showaverageanatomy',false,...
  'hfigs',hfigs_default,...
  'anatomydir',anatomydir,...
  'groupname',sprintf('Group selected %s',timestamp),...
  'int_manual',[]);
if numel(hfigs1) >= numel(hfigs_default),
  hfigs = hfigs1;
elseif isempty(hfigs1),
  hfigs = hfigs_default;
else
  hfigs = hfigs_default;
  hfigs(1:numel(hfigs1)) = hfigs1;
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

nlines = numel(line_names);

% create an html page with links to all experiments
filenamecurr = '';
if showbehavior || (showanatomy && ~isempty(int_manual)),
  filenamecurr = fullfile(tempdir,sprintf('selectedgroup_%s.html',timestamp));
  fid = fopen(filenamecurr,'w');
  if fid < 1,
    warning('Could not open temporary file %s for writing.\n');
  else
    
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
    
        
    if showbehavior,
      for i = 1:nlines,
        fprintf(fid,'<h1>%s</h1>\n',line_names{i});
        idx = find(strcmp({metadata.line_name},line_names{i}));
        if isempty(idx),
          fprintf(fid,'<p>No experiments found.</p>\n');
        else
          fprintf(fid,'<ul>\n');
          for j = idx(:)',
            [~,name] = fileparts(metadata(j).file_system_path);
            moviename = fullfile(metadata(j).file_system_path,...
              sprintf('ctrax_results_movie_%s.avi',name));
            plotsname = fullfile(metadata(j).file_system_path,'analysis_plots');
            if ~exist(moviename,'file') && ~exist(plotsname,'file'),
              continue;
            end
            fprintf(fid,'  <li>%s: ',name);
            if exist(moviename,'file'),
              fprintf(fid,'<a href="file://%s">Ctrax results movie</a>',moviename);
            end
            if exist(moviename,'file') && exist(plotsname,'dir'),
              fprintf(fid,', ');
            end
            if exist(plotsname,'dir'),
              fprintf(fid,'<a href="file://%s">Analysis plots</a>',plotsname);
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
          for j = idxanat(:)',            
            name = imdata(j).name;
            maxproj_file = imdata(j).maxproj_file_system_path;
            maxproj_ch2_file = regexprep(maxproj_file,'_total.jpg$','_ch2_total.jpg');
            reg_file = regexprep(maxproj_file,'_total.jpg$','.reg.local.jpg');
            translation_file = imdata(j).translation_file_path;
            if ~exist(maxproj_file,'file') && ~exist(translation_file,'file'),
              continue;
            end
            fprintf(fid,'  <li>%s: ',name);
            if exist(maxproj_file,'file'),
              fprintf(fid,'<a href="file://%s">Max projection image</a>',maxproj_file);
            end
            if exist(maxproj_ch2_file,'file') && ~strcmp(maxproj_file,maxproj_ch2_file),
              fprintf(fid,', <a href="file://%s">Max proj, channel 2</a>',maxproj_ch2_file);
            end
            if exist(reg_file,'file') && ~strcmp(maxproj_file,reg_file),
              fprintf(fid,', <a href="file://%s">Registered image</a>',reg_file);
            end
            if exist(translation_file,'file'),
              fprintf(fid,', <a href="file://%s">Translation video</a>',translation_file);
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
    
    fprintf(fid,'</body>\n</html>\n');
    fclose(fid);
    if ~exist(filenamecurr,'file'),
      warning('Could not open temporary file %s',filenamecurr);
    else
      % open this page
      web(filenamecurr,'-browser');
    end
  end
end

hax = [];
if showanatomy,
  
  if numel(hfigs) < 1,
    hfigs(1) = figure;
  end
  figure(hfigs(1));
  clf;
  nc = ceil(sqrt(nlines));
  nr = ceil((nlines)/nc);
  hax = createsubplots(nr,nc,[.01,.01],hfigs(1));
  if numel(hax) > nlines,
    delete(hax(nlines+1:end));
    hax = hax(1:nlines);
  end

  hwait = mywaitbar(0,'Reading max projection images');
  nlinesplotted = 0;
  for i = 1:nlines,
    filename = fullfile(anatomydir,sprintf('meanim_%s.png',line_names{i}));
    hwait = mywaitbar((i-1)/(nlines+1),hwait,sprintf('Reading %s...\n',line_names{i}));
    if exist(filename,'file'),
      im = imread(filename);
      image(im,'Parent',hax(i));
      axis(hax(i),'image');
      nlinesplotted = nlinesplotted+1;
    else
      fprintf('Line %s has no image\n',line_names{i});
    end
    title(hax(i),line_names{i},'Interpreter','none');
    set(hax(i),'XTick',[],'YTick',[]);
    drawnow;
  end

  if nlinesplotted > 0,
    impixelinfo(hfigs(1));
    linkaxes(hax);
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
  hax(end+1) = gca;
  imagesc(max(mu,[],3)','Parent',hax(nlines+1));
  axis(hax(nlines+1),'image');
  title(hax(nlines+1),'Mean image','Interpreter','none');
  colormap(hfigs(2),kjetsmooth(256));
  set(hax(nlines+1),'XTick',[],'YTick',[]);
  colorbar('peer',hax(nlines+1));
  impixelinfo(hfigs(2));
  linkaxes(hax);
  
  if ishandle(hwait),
    delete(hwait);
  end

end