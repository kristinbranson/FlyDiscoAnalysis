function ShowLineAnatomy(line_names,varargin)

persistent imdata;
persistent metadata;
persistent maskdata;
persistent cachedir;

timestamp = datestr(now,'yyyymmddTHHMMSS');

anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
defaultcachedir = '/nobackup/branson/AnatomyCacheData20131008';

defaultmetadatafile = 'CollectedPrimaryMetadata20130912.mat';
defaultimdatafile = 'ImageryData20130824.mat';
defaultmaskfile = 'FullBrainMask.mat';

[metadata,imdata,maskdata,anatomydir,cachedir,groupname] = myparse(varargin,...
  'metadata',metadata,...
  'imdata',imdata,...
  'maskdata',maskdata,...
  'anatomydir',anatomydir,...
  'cachedir',cachedir,...
  'groupname',sprintf('groupselected_%s',timestamp));

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
if isempty(metadata),

  [f,p] = uigetfile('*.mat','Select metadata mat file',defaultmetadatafile);
  if ~ischar(f),
    return;
  end
  metadatafile = fullfile(p,f);
  fprintf('Loading in metadata...\n');
  metadata = load(metadatafile,'metadata');
  metadata = metadata.metadata;
  
end

if isempty(imdata),
  [f,p] = uigetfile('*.mat','Select imdata mat file',defaultimdatafile);
  if ~ischar(f),
    return;
  end
  imagerydatafile = fullfile(p,f);
  fprintf('Loading in imagery data...\n');
  imdata = load(imagerydatafile,'imdata');
  imdata = imdata.imdata;
end  

if isempty(maskdata),
  [f,p] = uigetfile('*.mat','Select imdata mat file',defaultmaskfile);
  if ~ischar(f),
    return;
  end
  maskfile = fullfile(p,f);
  fprintf('Loading in mask data...');
  maskdata = load(maskfile);
  
  % symmetrize
  fprintf('Creating symmetric version of the mask...');
  mleft = regexp(maskdata.leg,'^(.*)_L$','once','tokens');
  isleft = ~cellfun(@isempty,mleft);
  idxleft = find(isleft);
  mleft = cellfun(@(x) x{1},mleft(isleft),'UniformOutput',false);
  mright = regexp(maskdata.leg,'^(.*)_R$','once','tokens');
  isright = ~cellfun(@isempty,mright);
  idxright = find(isright);
  mright = cellfun(@(x) x{1},mright(isright),'UniformOutput',false);
  isneither = ~isleft & ~isright;
  idxneither = find(isneither);
  maskdata.mask_symmetric = zeros(size(maskdata.mask),'uint8');
  ncompartments = numel(idxleft)+nnz(isneither);
  maskdata.leg_symmetric = cell(1,ncompartments);
  [~,idxleft2right] = ismember(mleft,mright);
  compi = 1;
  for i = 1:numel(idxneither),
    maskdata.mask_symmetric(maskdata.mask == idxneither(i)) = compi;
    maskdata.leg_symmetric{compi} = maskdata.leg{idxneither(i)};
    compi = compi + 1;
  end
  for i = 1:numel(idxleft),
    ileft = idxleft(i);
    iright = idxright(idxleft2right(i));
    maskdata.mask_symmetric(maskdata.mask == ileft | maskdata.mask == iright) = compi;
    maskdata.leg_symmetric{compi} = mleft{i};
    compi = compi + 1;
  end
  maskdata.mask_symmetric = permute(maskdata.mask_symmetric,[2,1,3]);
  
  fprintf('Computing spatial extent of each compartment...');
  % get the spatial extent of each compartment
  maskdata.xlims = nan(ncompartments,2);
  maskdata.ylims = nan(ncompartments,2);
  maskdata.zlims = nan(ncompartments,2);
  maskdata.npx = nan(1,ncompartments);
  for i = 1:ncompartments,
    is3 = any(any(maskdata.mask_symmetric==i,1),2);
    is2 = any(any(maskdata.mask_symmetric==i,1),3);
    is1 = any(any(maskdata.mask_symmetric==i,2),3);
    maskdata.xlims(i,:) = [find(is2,1),find(is2,1,'last')];
    maskdata.ylims(i,:) = [find(is1,1),find(is1,1,'last')];
    maskdata.zlims(i,:) = [find(is3,1),find(is3,1,'last')];
    maskdata.npx(i) = nnz(maskdata.mask_symmetric==i);
  end
  
end  
ncompartments = numel(maskdata.leg_symmetric);

%% what images are there per line?
imsperline = cell(1,nlines);
lineidx = cell(1,nlines);
for linei = 1:nlines,
  line_name = line_names{linei};  
  % find all the images of this line
  idx = find(strcmp({imdata.line},line_name));
  [~,order] = sort([imdata(idx).qi]);
  idx = idx(order);
  imsperline{linei} = {};  
  for jj = 1:numel(idx),
    j = idx(jj);
    imname = imdata(j).raw_file_system_path;
    if exist(imname,'file'),
      imsperline{linei}{end+1} = imname;
    else
      idx(jj) = 0;
    end
    lineidx{linei} = idx(idx~=0);
  end
end

%% what files exist and what files are missing?

% per-image, per-compartment files
perim_percomp_projfiles = cell(ncompartments+1,nlines);
perim_percomp_projmissing = cell(ncompartments+1,nlines);
for linei = 1:nlines,
  line_name = line_names{linei};
  for k = 1:ncompartments+1,
    perim_percomp_projmissing{linei,k} = true(1,numel(imsperline{line}));
    perim_percomp_projfiles{linei,k} = cell(1,numel(imsperline{line}));
    if k == ncompartments+1,
      compname = 'ALL';
    else
      compname = maskdata.leg_symmetric{k};
    end
    for j = 1:numel(imsperline{line}),
      [~,name] = fileparts(imsperline{line});
      name = regexp(name,'^(.*)\..*$','once','tokens');
      name = name{1};
      perim_percomp_projfiles{linei,k}{j} = fullfile(cachedir,line_name,sprintf('%s_%s_maxproj.png',compname,name));
      perim_percomp_projmissing{linei,k}(j) = ~exist(perim_percomp_projfiles{linei,k}{j},'file');
    end
  end
end

% per-line, per-compartment files
perline_percomp_projfiles = cell(ncompartments+1,nlines);
perline_percomp_projmissing = true(ncompartments+1,nlines);
for linei = 1:nlines,
  line_name = line_names{linei};
  for k = 1:ncompartments+1,
    if k == ncompartments+1,
      compname = 'ALL';
    else
      compname = maskdata.leg_symmetric{k};
    end
    perline_percomp_projfiles{linei,k} = fullfile(cachedir,line_name,sprintf('%s_%s_meanim_maxproj.png',compname,line_name));
    perline_percomp_projmissing(linei,k) = ~exist(perline_percomp_projfiles{linei,k},'file');
  end
end

% per-line text files
perline_txtfiles = cell(1,nlines);
perline_txtmissing = true(1,nlines);
for linei = 1:nlines,
  line_name = line_names{linei};
  perline_txtfiles{linei} = fullfile(cachedir,line_name,sprintf('%s_totalintensity_per_compartment.txt',line_name));
  perline_txtmissing(linei) = ~exist(perline_txtfiles{linei},'file');
end

% per-group, per-compartment files
pergroup_percomp_projfiles = cell(ncompartments+1,1);
pergroup_percomp_projmissing = true(ncompartments+1,1);
for k = 1:ncompartments+1,
  if k == ncompartments+1,
    compname = 'ALL';
  else
    compname = maskdata.leg_symmetric{k};
  end
  pergroup_percomp_projfiles{k} = fullfile(cachedir,groupname,sprintf('%s_%s_maxproj.png',compname,groupname));
  pergroup_percomp_projmissing(k) = ~exist(pergroup_percomp_projfiles{k},'file');
end
  
  
% compute the projections for each compartment in each image for each line
for linei = 1:nlines,
  line_name = line_names{linei};
  
  % find all the images of this line
  idx = find(strcmp({imdata.line},line_name));
  [~,order] = sort([imdata(idx).qi]);
  idx = idx(order);
  lineidx{linei} = idx;

  % create a directory to cache results if it doesn't exist already
  cachelinedir = fullfile(cachedir,line_name);
  if ~exist(cachelinedir,'dir'),
    mkdir(cachelinedir);
  end

  % for each image of this line
  for jj = 1:numel(idx),
    j = idx(jj);
    
    name = regexp(imdata(j).name,'.*/(.*)\..*','tokens','once');
    name = name{1};

    imname = imdata(j).raw_file_system_path;
    if isempty(imname) || ~exist(imname,'file'),
      continue;
    end
    
    fprintf('Computing per-compartment max-projections for %s...\n',name);
    
    filenames = cell(1,ncompartments+1);
    ismissing = false(1,ncompartments+1);
    for k = 1:ncompartments,
      filenames{k} = fullfile(cachelinedir,sprintf('%s_%s_maxproj.png',maskdata.leg_symmetric{k},name));
      ismissing(k) = ~exist(filenames{k},'file');
    end
    filenames{end} = fullfile(cachelinedir,sprintf('ALL_%s_maxproj.png',name));
    ismissing(end) = ~exist(filenames{end},'file');
    
    if any(ismissing),

      fprintf('Computing per-compartment max-projections for %s, line %d / %d, image %d / %d...\n',name,linei,nlines,jj,numel(idx));
      
      im = loadRaw2StackGreen(imname);
      
      idxmissing = find(ismissing);
      
      % compute per-compartment max projections
      for k = idxmissing(:)',
        if k == ncompartments + 1,
          immask = im;
          immask(maskdata.mask_symmetric==0) = 0;
        else        
          immask = im(maskdata.ylims(k,1):maskdata.ylims(k,2),...
            maskdata.xlims(k,1):maskdata.xlims(k,2),...
            maskdata.zlims(k,1):maskdata.zlims(k,2));
          immask(maskdata.mask_symmetric(maskdata.ylims(k,1):maskdata.ylims(k,2),...
            maskdata.xlims(k,1):maskdata.xlims(k,2),...
            maskdata.zlims(k,1):maskdata.zlims(k,2))~=k) = 0;
        end
        maxproj = max(immask,[],3)';
        maxproj = uint16(min(round(2^16*maxproj),2^16-1));
        imwrite(maxproj,filenames{k},'BitDepth',16);
      end      
      
    end
    
  end
    
end

% compute all-line mean image, per-compartment mean image projections

groupdir = fullfile(cachedir,groupname);
linenamefile = fullfile(groupdir,'line_names.txt');
meanimfile = fullfile(groupdir,'meanim.mat');
filenames = cell(1,ncompartments+1);
ismissing = false(1,ncompartments+1);
for k = 1:ncompartments,
  filenames{k} = fullfile(groupdir,sprintf('%s_%s_maxproj.png',maskdata.leg_symmetric{k},groupname));
  ismissing(k) = ~exist(filenames{k},'file');
end
filenames{end} = fullfile(groupdir,sprintf('ALL_%s_maxproj.png',groupname));
ismissing(end) = ~exist(filenames{end},'file');
    
isgroupdata = false;
if exist(groupdir,'dir') && ~any(ismissing) && exist(meanimfile,'file') && exist(linenamefile,'file'),
  
  line_names_read = importdata(linenamefile,',',0);
  if isempty(setdiff(line_names_read,line_names)) && isempty(setdiff(line_names,line_names_read)),
    isgroupdata = true;
  end
  
end

if ~isgroupdata,

[nr,nc,nz] = size(maskdata.mask_symmetric);
mu = zeros([nr,nc,nz],'single');
totalint_per_compartment = nan(nlines,ncompartments);

nread = 0;
for linei = 1:nlines,
  line_name = line_names{linei};
  
  imname = fullfile(anatomydir,sprintf('meanim_%s.png',line_name));
  if ~exist(imname,'file'),
    warning('File %s does not exist.',imname);
  end
  
  fprintf('Reading in mean image for line %s, %d/%d...\n',line_name,linei,nlines);
  im = load(imname,'meanim');
  im = im.meanim;
  mu = mu + im;
  nread = nread + 1;

  cachelinedir = fullfile(cachedir,line_name);
  if ~exist(cachelinedir,'dir'),
    mkdir(cachelinedir);
  end

  filenames = cell(1,ncompartments);
  ismissing = false(1,ncompartments);
  for k = 1:ncompartments,
    filenames{k} = fullfile(cachelinedir,sprintf('%s_%s_meanim_maxproj.png',maskdata.leg_symmetric{k},line_name));
    ismissing(k) = ~exist(filenames{k},'file');
  end
  
  txtfile = fullfile(cachelinedir,sprintf('%s_totalintensity_per_compartment.txt',line_name));
  txtfilemissing = ~exist(txtfile,'file');
  
  if any(ismissing) || txtfilemissing,
    
    fprintf('Computing per-compartment max-projections for mean image for %s, line %d / %d...\n',line_name,linei,nlines);
    idxmissing = find(ismissing);
    if txtfilemissing,
      idxmissing = 1:ncompartments;
    end
    
    % compute per-compartment max projections
    for k = idxmissing(:)',
      immask = im(maskdata.ylims(k,1):maskdata.ylims(k,2),...
        maskdata.xlims(k,1):maskdata.xlims(k,2),...
        maskdata.zlims(k,1):maskdata.zlims(k,2));
      immask(maskdata.mask_symmetric(maskdata.ylims(k,1):maskdata.ylims(k,2),...
        maskdata.xlims(k,1):maskdata.xlims(k,2),...
        maskdata.zlims(k,1):maskdata.zlims(k,2))~=k) = 0;
      if ismissing(k),
        maxproj = max(immask,[],3)';
        maxproj = uint16(min(round(2^16*maxproj),2^16-1));
        imwrite(maxproj,filenames{k},'BitDepth',16);
      end
      totalint_per_compartment(linei,k) = sum(single(immask(:))) / npx;
    end

    if txtfilemissing,
      fid = fopen(txtfile,'w');
      for k = 1:ncompartments,
        fprintf(fid,'%s,%f\n',maskdata.leg_symmetric{k},totalint_per_compartment(linei,k));
      end
      fclose(fid);
    end
    
  end
  
  if ~txtfilemissing,
    tmp = importdata(txtfile,',',0);
    [~,compidx] = ismember(tmp(:,1),maskdata.leg_symmetric);
    totalint_per_compartment(linei,compidx) = cell2mat(tmp(:,2));
  end
  
end

mu = mu / nread;

% compute per-compartment mean images
for k = 1:ncompartments,
  
end