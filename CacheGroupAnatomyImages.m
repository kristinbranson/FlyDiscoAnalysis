function CacheGroupAnatomyImages(line_names,varargin)

persistent maskdata;
persistent cachedir;
persistent alllines_meanim;

anatomydir = '/nobackup/branson/AverageAnatomyData20130618';
defaultcachedir = '/nobackup/branson/AnatomyCacheData20131008';
defaultalllinematfile = 'alllines_anatomy_20130617/meanim.mat';

defaultmaskfile = 'FullBrainMaskSymmetric.mat';
timestamp = datestr(now,'yyyymmddTHHMMSS');

[maskdata,anatomydir,cachedir,maskfile,groupname,alllines_meanim,alllinematfile] = ...
  myparse(varargin,...
  'maskdata',maskdata,...
  'anatomydir',anatomydir,...
  'cachedir',cachedir,...
  'maskfile','',...
  'groupname',sprintf('selected_%s',timestamp),...
  'alllines_meanim',alllines_meanim,...
  'alllinematfile','');

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

if isempty(maskdata),
  if isempty(maskfile),
    [f,p] = uigetfile('*.mat','Select maskdata mat file',defaultmaskfile);
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
  
if isempty(alllines_meanim),
  
  if isempty(alllinematfile),
    [f,p] = uigetfile('*.mat','Select all line mean image mat file',defaultalllinematfile);
    if ~ischar(f),
      return;
    end
    alllinematfile = fullfile(p,f);
  end
  fprintf('Loading in all line mean image...');
  alllines_meanim = load(alllinematfile);
  alllines_meanim = alllines_meanim.meanim;

end

%% compute the mean image
nread = 0;
nlines = numel(line_names);
for linei = 1:nlines,
  
  line_name = line_names{linei};
  imname = fullfile(anatomydir,sprintf('meanim_%s.mat',line_name));
  fprintf('Reading in mean image for line %s, %d/%d...\n',line_name,linei,nlines);
  if ~exist(imname,'file'),
    warning('File %s does not exist',imname);
    continue;
  end
  load(imname,'meanim');
  
  if linei == 1,
    mu = meanim;
  else
    mu = mu + meanim;
  end
  
  nread = nread + 1;
  
end

mu = mu / nread;

%% compute the difference from the mean

mudiff = mu - alllines_meanim;

%% create the directory

cachegroupdir = fullfile(cachedir,groupname);
if ~exist(cachegroupdir,'dir'),
  mkdir(cachegroupdir);
end

%% save the mean image to a mat file

fprintf('Saving mean image to mat file...\n');
filename = fullfile(cachegroupdir,sprintf('%s_meanim.mat',groupname));
save('-v7.3',filename,'mu');

%% save the mean image projection images

totalint_per_compartment = nan(1,ncompartments);  

filenames = cell(1,ncompartments+1);
for k = 1:ncompartments,
  filenames{k} = fullfile(cachegroupdir,sprintf('%s_%s_meanim_maxproj.png',maskdata.leg_symmetric{k},groupname));
end
filenames{end} = fullfile(cachegroupdir,sprintf('ALL_%s_meanim_maxproj.png',groupname));
txtfile = fullfile(cachegroupdir,sprintf('%s_totalintensity_per_compartment.txt',groupname));

% compute per-compartment max projections
maxv = nan(1,ncompartments+1);
for k = 1:ncompartments+1,
  if k == ncompartments + 1,
    immask = mu;
    immask(maskdata.mask_symmetric==0) = 0;
  else
    immask = mu(maskdata.ylims(k,1):maskdata.ylims(k,2),...
      maskdata.xlims(k,1):maskdata.xlims(k,2),...
      maskdata.zlims(k,1):maskdata.zlims(k,2));
    immask(maskdata.mask_symmetric(maskdata.ylims(k,1):maskdata.ylims(k,2),...
      maskdata.xlims(k,1):maskdata.xlims(k,2),...
      maskdata.zlims(k,1):maskdata.zlims(k,2))~=k) = 0;
    totalint_per_compartment(k) = sum(single(immask(:)));
  end
  maxv(k) = max(immask(:));
  maxproj = max(immask,[],3)';
  maxproj = uint16(min(round(2^16*maxproj/maxv(k)),2^16-1));
  imwrite(maxproj,filenames{k},'BitDepth',16);
end

fid = fopen(txtfile,'w');
for k = 1:ncompartments,
  fprintf(fid,'%s,%f,%f\n',maskdata.leg_symmetric{k},totalint_per_compartment(k),maxv(k));
end
fprintf(fid,'ALL,0,%f\n',maxv(ncompartments+1));
fclose(fid);

%% save the mean image difference projection images

totalint_per_compartment = nan(1,ncompartments);  

filenames = cell(1,ncompartments+1);
for k = 1:ncompartments,
  filenames{k} = fullfile(cachegroupdir,sprintf('%s_%s_meanimdiff_maxproj.png',maskdata.leg_symmetric{k},groupname));
end
filenames{end} = fullfile(cachegroupdir,sprintf('ALL_%s_meanimdiff_maxproj.png',groupname));
txtfile = fullfile(cachegroupdir,sprintf('%s_diff_totalintensity_per_compartment.txt',groupname));

% compute per-compartment max projections
maxv = nan(1,ncompartments+1);
minv = nan(1,ncompartments+1);
for k = 1:ncompartments+1,
  if k == ncompartments + 1,
    immask = mudiff;
    immask(maskdata.mask_symmetric==0) = 0;
  else
    immask = mudiff(maskdata.ylims(k,1):maskdata.ylims(k,2),...
      maskdata.xlims(k,1):maskdata.xlims(k,2),...
      maskdata.zlims(k,1):maskdata.zlims(k,2));
    immask(maskdata.mask_symmetric(maskdata.ylims(k,1):maskdata.ylims(k,2),...
      maskdata.xlims(k,1):maskdata.xlims(k,2),...
      maskdata.zlims(k,1):maskdata.zlims(k,2))~=k) = 0;
    totalint_per_compartment(k) = sum(single(immask(:)));
  end
  maxproj = max(immask,[],3)';
  maxv(k) = max(maxproj(:));
  minv(k) = min(maxproj(:));
  maxproj = uint16(min(round(2^16*(maxproj-minv(k))/(maxv(k)-minv(k))),2^16-1));
  imwrite(maxproj,filenames{k},'BitDepth',16);
end

fid = fopen(txtfile,'w');
for k = 1:ncompartments,
  fprintf(fid,'%s,%f,%f,%f\n',maskdata.leg_symmetric{k},totalint_per_compartment(k),maxv(k),minv(k));
end
fprintf(fid,'ALL,0,%f,%f\n',maxv(ncompartments+1),minv(ncompartments+1));
fclose(fid);

%% save a text file with the line names for this group

filename = fullfile(cachegroupdir,sprintf('%s_line_names.txt',groupname));
fid = fopen(filename,'w');
for linei = 1:nlines,
  fprintf(fid,'%s\n',line_names{linei});
end
fclose(fid);
