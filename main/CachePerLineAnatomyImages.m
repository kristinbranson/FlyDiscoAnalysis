function CachePerLineAnatomyImages(line_name,varargin)

persistent imdata;
persistent maskdata;
persistent cachedir;

anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
defaultcachedir = '/nobackup/branson/AnatomyCacheData20131008';

defaultimdatafile = 'ImageryData20130824.mat';
defaultmaskfile = 'FullBrainMaskSymmetric.mat';

[imdata,maskdata,anatomydir,cachedir,maxqi,imdatafile,maskfile] = myparse(varargin,...
  'imdata',imdata,...
  'maskdata',maskdata,...
  'anatomydir',anatomydir,...
  'cachedir',cachedir,...
  'maxqi',.3,...
  'imdatafile','',...
  'maskfile','');

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
  
  % symmetrize
  maskdata = SymmetrizeMaskData(maskdata);
  
end  
ncompartments = numel(maskdata.leg_symmetric);  
  
% compute the projections for each compartment in each image of the line

% find all the images of this line
idx = find(strcmp({imdata.line},line_name));
idx([imdata(idx).qi]>maxqi) = [];
[~,order] = sort([imdata(idx).qi]);
idx = idx(order);

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

% compute per-compartment mean image projections
totalint_per_compartment = nan(1,ncompartments);
imname = fullfile(anatomydir,sprintf('meanim_%s.mat',line_name));
if ~exist(imname,'file'),
  warning('File %s does not exist.',imname);
else
  fprintf('Reading in mean image for line %s...\n',line_name);
  im = load(imname,'meanim');
  im = im.meanim;
  
  filenames = cell(1,ncompartments+1);
  ismissing = false(1,ncompartments+1);
  for k = 1:ncompartments,
    filenames{k} = fullfile(cachelinedir,sprintf('%s_%s_meanim_maxproj.png',maskdata.leg_symmetric{k},line_name));
    ismissing(k) = ~exist(filenames{k},'file');
  end
  filenames{end} = fullfile(cachelinedir,sprintf('ALL_%s_meanim_maxproj.png',line_name));
  ismissing(end) = ~exist(filenames{end},'file');
  
  txtfile = fullfile(cachelinedir,sprintf('%s_totalintensity_per_compartment.txt',line_name));
  txtfilemissing = ~exist(txtfile,'file');
  
  if any(ismissing) || txtfilemissing,
    
    fprintf('Computing per-compartment max-projections for mean image for %s...\n',line_name);
    idxmissing = find(ismissing);
    if txtfilemissing,
      idxmissing = 1:ncompartments+1;
    end
    
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
        totalint_per_compartment(k) = sum(single(immask(:)));
      end
      if ismissing(k),
        maxproj = max(immask,[],3)';
        maxproj = uint16(min(round(2^16*maxproj),2^16-1));
        imwrite(maxproj,filenames{k},'BitDepth',16);
      end
    end

    if txtfilemissing,
      fid = fopen(txtfile,'w');
      for k = 1:ncompartments,
        fprintf(fid,'%s,%f\n',maskdata.leg_symmetric{k},totalint_per_compartment(k));
      end
      fclose(fid);
    end
    
  end
    
end