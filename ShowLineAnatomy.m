function ShowLineAnatomy(line_names,varargin)

persistent imdata;
persistent metadata;
persistent maskdata;
persistent cachedir;

tmpstruct = struct;

timestamp = datestr(now,'yyyymmddTHHMMSS');

anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
defaultcachedir = '/nobackup/branson/AnatomyCacheData20131008';

defaultmetadatafile = 'CollectedPrimaryMetadata20130912.mat';
defaultimdatafile = 'ImageryData20130824.mat';
defaultmaskfile = 'FullBrainMaskSymmetric.mat';

controlwidth = 150;
controlrightborder = 20;

[metadata,imdata,maskdata,anatomydir,cachedir,groupname,...
  metadatafile,imdatafile,maskfile,hfigs,figpos,maxqi] = myparse(varargin,...
  'metadata',metadata,...
  'imdata',imdata,...
  'maskdata',maskdata,...
  'anatomydir',anatomydir,...
  'cachedir',cachedir,...
  'groupname',sprintf('groupselected_%s',timestamp),...
  'metadatafile','',...
  'imdatafile','',...
  'maskfile','',...
  'hfigs',[],...
  'figpos',get(0,'ScreenSize'),...
  'maxqi',.3);

nlines = numel(line_names);

%% persistent data

if isempty(cachedir),
  p = uigetdir(defaultcachedir,'Select directory to cache results to');
  if ~ischar(p),
    return;
  end
  cachedir = p;
end

if ~exist(cachedir,'dir'),
  mkdir(cachedir);
end

% load in missing data
% if isempty(metadata),
% 
%   if isempty(metadatafile),
%     [f,p] = uigetfile('*.mat','Select metadata mat file',defaultmetadatafile);
%     if ~ischar(f),
%       return;
%     end
%     metadatafile = fullfile(p,f);
%   end
%   fprintf('Loading in metadata...\n');
%   metadata = load(metadatafile,'metadata');
%   metadata = metadata.metadata;
%   
% end

if isempty(imdata),
  if isempty(imdatafile),
    [f,p] = uigetfile('*.mat','Select imdata mat file',defaultimdatafile);
    if ~ischar(f),
      return;
    end
    imdatafile = fullfile(p,f);
  end
  fprintf('Loading in imagery data...\n');
  imdata = load(imdatafile,'imdata');
  imdata = imdata.imdata;
end  

if isempty(maskdata),
  if isempty(maskfile),
    [f,p] = uigetfile('*.mat','Select imdata mat file',defaultmaskfile);
    if ~ischar(f),
      return;
    end
    maskfile = fullfile(p,f);
  end
  fprintf('Loading in mask data...');
  maskdata = load(maskfile);
  
  maskdata = SymmetrizeMaskData(maskdata);
end  
ncompartments = numel(maskdata.leg_symmetric);
colors = jet(ncompartments);

%% what images are cached per-line

imnames_perline = cell(1,nlines);
nims_perline = nan(1,nlines);
maxv_perline_perim = cell(1,nlines);
for linei = 1:nlines,

  line_name = line_names{linei};
  cachelinedir = fullfile(cachedir,line_name);
  if ~exist(cachelinedir,'dir'),
    warning('Cache directory for line %s, %s, does not exist',line_name,cachelinedir);
    continue;
  end

  idx = find(strcmp({imdata.line},line_name));
  idx([imdata(idx).qi]>maxqi) = [];
  [~,order] = sort([imdata(idx).qi]);
  idx = idx(order);
  
  imnames_perline{linei} = {};
  maxv_perline_perim{linei} = [];
  for jj = 1:numel(idx),
    j = idx(jj);
    name = regexp(imdata(j).name,'.*/(.*)\..*','tokens','once');
    name = name{1};
    filename = fullfile(cachelinedir,sprintf('ALL_%s_maxproj.png',name));
    if exist(filename,'file'),
      imnames_perline{linei}{end+1} = name;
    end
    im = single(imread(filename))/(2^16-1);
    maxv_perline_perim{linei}(end+1) = max(im(:));
  end
  nims_perline(linei) = numel(imnames_perline{linei});
end

%% create the figure

if isempty(hfigs),
  hfigs(1) = figure;
else
  figure(hfigs(1));
  clf;
  set(hfigs(1),'Units','pixels');
end
set(hfigs(1),'Position',figpos);

%% create the axes for each line

linei = 1;
line_name = line_names{linei};
cachelinedir = fullfile(cachedir,line_name);
filename = fullfile(cachelinedir,sprintf('ALL_%s_meanim_maxproj.png',line_name));
im = imread(filename);
[nr,nc] = size(im);
compartment_names = maskdata.leg_symmetric;
ncompartments = numel(compartment_names);

maxscale = 0;
for ncax = 3:nlines,
  nrax = ceil(nlines/ncax);
  scale = min((figpos(3)-controlwidth+controlrightborder)/(nc*ncax),figpos(4)/((nr+1)*nrax));
  if scale > maxscale,
    maxscale = scale;
    bestnrax = nrax;
  end
end
nrax = bestnrax;
ncax = ceil(nlines/nrax);

hax = nan(nrax,ncax);
w = (figpos(3)-controlwidth)/figpos(3)/ncax;
h = 1/nrax;
f = .02;
for r = 1:nrax,
  offr = r*h + h*f/2;
  for c = 1:ncax,
    offc = (c-1)*w + w*f/2;
    hax(r,c) = axes('Position',[offc,offr,w*(1-f),h*(1-f)]); %#ok<LAXES>
  end
end
hax = flipud(hax);

if nrax*ncax > nlines,
  delete(hax(nlines+1:end));
end

%% show the per-line images

shortlinenames = regexprep(line_names,'^GMR_','R');
shortlinenames = regexprep(shortlinenames,'_01','');
shortlinenames = regexprep(shortlinenames,'_AE','');

colormap(kjetsmooth(256));

him = nan(1,nlines);
for linei = 1:nlines,
  
  line_name = line_names{linei};
  cachelinedir = fullfile(cachedir,line_name);
  if ~exist(cachelinedir,'dir'),
    continue;
  end
  
  filename = fullfile(cachelinedir,sprintf('ALL_%s_meanim_maxproj.png',line_name));
  if ~exist(filename,'file'),
    continue;
  end
  im = single(imread(filename))/(2^16-1);
  him(linei) = imagesc(im,'Parent',hax(linei),[0,1]);
  axis(hax(linei),'image','off');
  
  text(0,0,shortlinenames{linei},'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','w','Parent',hax(linei));
  
  drawnow;
end

%% create the listbox

offx = ncax*w + f/2*w;
offy = f+.25;
pos = [offx,offy,1-offx-2*f*w,.75-2*f];

hcompartments = uicontrol('Style','Listbox','Units','normalized','Position',pos,'Max',10);
compartment_list = [compartment_names,{'All'}];
set(hcompartments,'String',compartment_list,'Value',ncompartments+1);
compartments_selected = {'All'};
is_compartment_selected = true(1,ncompartments);

offy = offy - 40/figpos(3);
pos = [offx,offy,1-offx-2*f*w,20/figpos(3)];
show_list = {'Mean image','Single image'};
hshow = uicontrol('Style','popupmenu','Units','normalized','Position',pos,'String',show_list,'Value',1);
show_mode = 'Mean image';
imidx_perline = ones(1,nlines);

offy = offy - 80/figpos(3);
pos = [offx,offy,1-offx-2*f*w,40/figpos(3)];
hupdate = uicontrol('Style','pushbutton','Units','normalized','Position',pos,'String','Update','Callback',@UpdateLineAnatomy);

% put a button on each image to cycle through lines
hnextim = nan(1,nlines);
for linei = 1:nlines,
  pos = get(hax(linei),'Position');
  pos(3:4) = [30/figpos(3),30/figpos(4)];
  hnextim(linei) = uicontrol('Style','pushbutton','Units','normalized','Position',pos,'String','>','Callback',{@NextImage,linei});
end

offy = offy - 80/figpos(3);
pos = [offx,offy,1-offx-2*f*w,40/figpos(3)];
hshowcomps = uicontrol('Style','checkbox','Units','normalized','Position',pos,'String','Show compartments','Callback',@ToggleShowCompartments);


%% look for group data

cachegroupdir = fullfile(cachedir,groupname);
isgroupdata = false;
if exist(cachegroupdir,'dir'),
  line_name_file = fullfile(cachegroupdir,sprintf('%s_line_names.txt',groupname));
  if exist(line_name_file,'file'),
    line_names_read = importdata(line_name_file);
    if isempty(setdiff(line_names_read,line_names)) && isempty(setdiff(line_names,line_names_read)),
      isgroupdata = true;
    end
  end
end

if ~isgroupdata,
  drawnow;
  fprintf('Computing group mean images...\n');
  CacheGroupAnatomyImages(line_names,'maskdata',maskdata,'anatomydir',anatomydir,'cachedir',cachedir,...
    'groupname',groupname);
end

% read in the compartment text files
txtfile = fullfile(cachegroupdir,sprintf('%s_totalintensity_per_compartment.txt',groupname));
totalint_per_compartment_mean = nan(1,ncompartments);
maxv_per_compartment_mean = nan(1,ncompartments+1);
fid = fopen(txtfile,'r');
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  s = regexp(s,',','split');
  if strcmp(s{1},'ALL'),
    maxv_per_compartment_mean(ncompartments+1) = str2double(s{3});
  else
    k = find(strcmp(compartment_names,s{1}));
    totalint_per_compartment_mean(k) = str2double(s{2});
    maxv_per_compartment_mean(k) = str2double(s{3});    
  end
end
fclose(fid);

% read in the compartment text file for the meandiff
txtfile = fullfile(cachegroupdir,sprintf('%s_diff_totalintensity_per_compartment.txt',groupname));
totalint_per_compartment_meandiff = nan(1,ncompartments);
maxv_per_compartment_meandiff = nan(1,ncompartments+1);
minv_per_compartment_meandiff = nan(1,ncompartments+1);
fid = fopen(txtfile,'r');
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  s = regexp(s,',','split');
  if strcmp(s{1},'ALL'),
    maxv_per_compartment_meandiff(ncompartments+1) = str2double(s{3});
    minv_per_compartment_meandiff(ncompartments+1) = str2double(s{4});
  else
    k = find(strcmp(compartment_names,s{1}));
    totalint_per_compartment_meandiff(k) = str2double(s{2});
    maxv_per_compartment_meandiff(k) = str2double(s{3});    
    minv_per_compartment_meandiff(k) = str2double(s{4});
  end
end
fclose(fid);

%% create group data axes

hax_group = nan(1,3);
offr = h*f/2;
w = (figpos(3)-controlwidth)/figpos(3)/3;
for c = 1:3,
  offc = (c-1)*w + w*f/2;
  hax_group(c) = axes('Position',[offc,offr,w*(1-f),h*(1-f)]); %#ok<LAXES>
end

%% show mean images

filename = fullfile(cachegroupdir,sprintf('ALL_%s_meanim_maxproj.png',groupname));
im = single(imread(filename))/(2^16-1) * maxv_per_compartment_mean(ncompartments+1);
him_group = nan(1,2);
him_group(1) = imagesc(im,'Parent',hax_group(1),[0,maxv_per_compartment_mean(ncompartments+1)]);
text(0,0,'mean','Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','w','Parent',hax_group(1));
axis(hax_group(1),'image','off');
colorbar('peer',hax_group(1));

filename = fullfile(cachegroupdir,sprintf('ALL_%s_meanimdiff_maxproj.png',groupname));
im = single(imread(filename))/(2^16-1);
minv = minv_per_compartment_meandiff(ncompartments+1);
maxv = maxv_per_compartment_meandiff(ncompartments+1);
im = im*(maxv-minv)+minv;
him_group(2) = imagesc(im,'Parent',hax_group(2),[0,maxv]);
text(0,0,'meandiff','Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','w','Parent',hax_group(2));
axis(hax_group(2),'image','off');
colorbar('peer',hax_group(2));

linkaxes([hax(ishandle(hax));hax_group(1:2)']);

impixelinfo;

drawnow;

%% plot the compartment outlines
hlines = nan(nlines+2,ncompartments);
for linei = 1:nlines,
  for k = 1:ncompartments,
    hold(hax(linei),'on');
    hlines(linei,k) = plot(hax(linei),maskdata.boundaries(k).y,maskdata.boundaries(k).x,'-','Color',colors(k,:));
  end
end
for i = 1:2,
  for k = 1:ncompartments,
    hold(hax_group(i),'on');
    hlines(nlines+i,k) = plot(hax_group(i),maskdata.boundaries(k).y,maskdata.boundaries(k).x,'-','Color',colors(k,:));
  end
end  

set(hlines,'Visible','off');

%% show compartment average intensities

% totalint_per_compartment = nan(ncompartments,nlines+2);
% for linei = 1:nlines,
%   line_name = line_names{linei};
%   cachelinedir = fullfile(cachedir,line_name);
%   txtfile = fullfile(cachelinedir,sprintf('%s_totalintensity_per_compartment.txt',line_name));
%   
%   fid = fopen(txtfile,'r');
%   while true,
%     s = fgetl(fid);
%     if ~ischar(s),
%       break;
%     end
%     s = strtrim(s);
%     s = regexp(s,',','split');
%     if strcmp(s{1},'ALL'),
%       continue;
%     end
%     k = strcmp(compartment_names,s{1});
%     totalint_per_compartment(k,linei) = str2double(s{2});
%   end
%   fclose(fid);
% end
% 
% totalint_per_compartment(:,nlines+1) = totalint_per_compartment_mean';
% totalint_per_compartment(:,nlines+2) = totalint_per_compartment_meandiff';
% meanint_per_compartment = bsxfun(@rdivide,totalint_per_compartment,maskdata.npx(:));
meanint_per_compartment = totalint_per_compartment_mean ./ maskdata.npx;

off = 1;
for k = 1:ncompartments,
  x = [off+.05,off+.95];
  y = [0,meanint_per_compartment(k)];
  color = colors(k,:);
  patch(x([1,1,2,2]),y([1,2,2,1]),color,'Parent',hax_group(3),'EdgeColor','w');
  off = off + 1;
end
axis(hax_group(3),[0,off,min(meanint_per_compartment(:)),max(meanint_per_compartment(:))]);
set(hax_group(3),'XTick',1.5:ncompartments+.5,'XTickLabel',compartment_names,'Color','k');
pos = get(hax_group(3),'Position');
tmpstruct.addon = pos(4)*.2;
pos(2) = pos(2) + tmpstruct.addon;
pos(4) = pos(4) - tmpstruct.addon;
set(hax_group(3),'Position',pos);
rotateticklabel(hax_group(3));
ylabel(hax_group(3),'Mean mean');

%% callback

  function UpdateLineAnatomy(varargin)

    old_compartments_selected = compartments_selected;
    old_show_mode = show_mode;
    
    % update compartments selected
    v = get(hcompartments,'Value');
    compartments_selected_curr = compartment_list(v);
    if ismember('All',compartments_selected_curr),
      set(hcompartments,'Value',ncompartments+1);
      compartments_selected = {'All'};
    elseif isempty(compartments_selected_curr),
      [~,v] = ismember(compartments_selected,compartment_list);
      set(hcompartments,'Value',v);
      warning('At least one compartment must be selected');
      return;
    else
      compartments_selected = compartments_selected_curr;
    end
    
    % update plot mode
    v = get(hshow,'Value');
    show_mode = show_list{v};

    ischange = ~isempty(setdiff(old_compartments_selected,compartments_selected)) || ...
      ~isempty(setdiff(compartments_selected,old_compartments_selected)) || ...
      ~strcmp(show_mode,old_show_mode);
    
    if ~ischange,
      return;
    end
    
    if ismember('All',compartments_selected),
      is_compartment_selected = true(1,ncompartments);
    else
      is_compartment_selected = ismember(maskdata.leg_symmetric,compartments_selected);
    end
    
    % update the images
    for linei1 = 1:nlines,
      line_name = line_names{linei1};
      cachelinedir = fullfile(cachedir,line_name);
      maxim = zeros([nr,nc],'single');
      for k = 1:numel(compartments_selected),
        compartment = upper(compartments_selected{k});
        if strcmp(show_mode,'Mean image'),
          filename = fullfile(cachelinedir,sprintf('%s_%s_meanim_maxproj.png',compartment,line_name));
        else
          name = imnames_perline{linei1}{imidx_perline(linei1)};
          filename = fullfile(cachelinedir,sprintf('%s_%s_maxproj.png',compartment,name));
        end
        if ~exist(filename,'file'),
          warning('File %s does not exist, not changing image for this line.',filename);
          continue;
        end
        im = single(imread(filename))/(2^16-1);
        if ~strcmp(compartment,'ALL'),
          kk = find(strcmp(compartment,maskdata.leg_symmetric));
          xlim = maskdata.xlims(kk,:);
          ylim = maskdata.ylims(kk,:);
          %zlim = maskdata.zlims(kk,:);
          % note that we switch x and y because the transpose is taken
          % while saving
          maxim(xlim(1):xlim(2),ylim(1):ylim(2)) = ...
            max(maxim(xlim(1):xlim(2),ylim(1):ylim(2)),im);
        else
          maxim = max(maxim,im);
        end
      end
      set(him(linei1),'CData',maxim);
      if strcmp(show_mode,'Mean image'),
        set(hax(linei1),'CLim',[0,1]);
      else
        set(hax(linei1),'CLim',[0,maxv_perline_perim{linei1}(imidx_perline(linei1))]);
      end
    end
    
    % update the mean images
    maxim = zeros([nr,nc],'single');
    for k = 1:numel(compartments_selected),
      compartment = upper(compartments_selected{k});
      
      filename = fullfile(cachegroupdir,sprintf('%s_%s_meanim_maxproj.png',compartment,groupname));
      im = single(imread(filename))/(2^16-1);
      
      if ~strcmp(compartment,'ALL'),
        kk = find(strcmp(compartment,maskdata.leg_symmetric));
        maxv = maxv_per_compartment_mean(kk);
        im = im*maxv;
        xlim = maskdata.xlims(kk,:);
        ylim = maskdata.ylims(kk,:);
        % note that we switch x and y because the transpose is taken
        % while saving
        maxim(xlim(1):xlim(2),ylim(1):ylim(2)) = ...
          max(maxim(xlim(1):xlim(2),ylim(1):ylim(2)),im);
      else
        maxv = maxv_per_compartment_mean(ncompartments+1);
        im = im*maxv;
        maxim = max(maxim,im);
      end
    end
    set(him_group(1),'CData',maxim);
    set(hax_group(1),'CLim',[0,max(maxim(:))]);
    
    maxim = zeros([nr,nc],'single');
    for k = 1:numel(compartments_selected),
      compartment = upper(compartments_selected{k});
      
      filename = fullfile(cachegroupdir,sprintf('%s_%s_meanimdiff_maxproj.png',compartment,groupname));
      im = single(imread(filename))/(2^16-1);
      
      if ~strcmp(compartment,'ALL'),
        kk = find(strcmp(compartment,maskdata.leg_symmetric));
        maxv = maxv_per_compartment_meandiff(kk);
        minv = minv_per_compartment_meandiff(kk);
        im = im*(maxv-minv)+minv;
        xlim = maskdata.xlims(kk,:);
        ylim = maskdata.ylims(kk,:);
        % note that we switch x and y because the transpose is taken
        % while saving
        maxim(xlim(1):xlim(2),ylim(1):ylim(2)) = ...
          max(maxim(xlim(1):xlim(2),ylim(1):ylim(2)),im);
      else
        maxv = maxv_per_compartment_mean(ncompartments+1);
        minv = minv_per_compartment_meandiff(ncompartments+1);
        im = im*(maxv-minv)+minv;
        maxim = max(maxim,im);
      end
    end
    set(him_group(2),'CData',maxim);
    set(hax_group(2),'CLim',[0,max(maxim(:))]);
    
    % update whether the compartment outlines are shown
    ToggleShowCompartments();

      
% him_group = nan(1,2);
% him_group(1) = imagesc(im,'Parent',hax_group(1),[0,maxv_per_compartment_mean(ncompartments+1)]);
% text(0,0,'mean','Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','w','Parent',hax_group(1));
% 
% axis(hax_group(1),'image','off');
% 
% filename = fullfile(cachegroupdir,sprintf('ALL_%s_meanimdiff_maxproj.png',groupname));
% im = single(imread(filename))/(2^16-1);
% minv = minv_per_compartment_meandiff(ncompartments+1);
% maxv = maxv_per_compartment_meandiff(ncompartments+1);
% im = im*(maxv-minv)+minv;
% him_group(2) = imagesc(im,'Parent',hax_group(2),[0,maxv]);
% text(0,0,'meandiff','Interpreter','none','HorizontalAlignment','left','VerticalAlignment','top','Color','w','Parent',hax_group(2));
% axis(hax_group(2),'image','off');
% 
    
    
  end

  function NextImage(~,~,linei1)
    
    if strcmp(show_mode,'Mean image'),
      warning('Cannot change image while in "Mean image" mode');
      return;
    end
    imidx_perline(linei1) = imidx_perline(linei1) + 1;
    if imidx_perline(linei1) > nims_perline(linei1),
      imidx_perline(linei1) = 1;
    end
    
    line_name = line_names{linei1};
    cachelinedir = fullfile(cachedir,line_name);
    maxim = zeros([nr,nc],'single');
    for k = 1:numel(compartments_selected),
      compartment = upper(compartments_selected{k});
      name = imnames_perline{linei1}{imidx_perline(linei1)};
      filename = fullfile(cachelinedir,sprintf('%s_%s_maxproj.png',compartment,name));
      if ~exist(filename,'file'),
        warning('File %s does not exist, not changing image for this line.',filename);
        continue;
      end
      im = single(imread(filename))/(2^16-1);
      if ~strcmp(compartment,'ALL'),
        kk = find(strcmp(compartment,maskdata.leg_symmetric));
        xlim = maskdata.xlims(kk,:);
        ylim = maskdata.ylims(kk,:);
        %zlim = maskdata.zlims(kk,:);
        % note that we switch x and y because the transpose is taken
        % while saving
        maxim(xlim(1):xlim(2),ylim(1):ylim(2)) = ...
          max(maxim(xlim(1):xlim(2),ylim(1):ylim(2)),im);
      else
        maxim = max(maxim,im);
      end
    end
    set(him(linei1),'CData',maxim);
    set(hax(linei1),'CLim',[0,maxv_perline_perim{linei1}(imidx_perline(linei1))]);
  end

  function ToggleShowCompartments(varargin)
    
    v = get(hshowcomps,'Value');
    if v,
      set(hlines(:,is_compartment_selected),'Visible','on');
      set(hlines(:,~is_compartment_selected),'Visible','off');
    else
      set(hlines,'Visible','off');
    end
    
  end

end