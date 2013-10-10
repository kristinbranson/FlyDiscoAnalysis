
imdatafile = 'ImageryData20130824.mat';
maskfile = 'FullBrainMaskSymmetric.mat';
cachedir = '/nobackup/branson/AnatomyCacheData20131008';
alllinematfile = 'alllines_anatomy_20130617/meanim.mat';

line_name = 'GMR_44D11_AE_01';
figpos = [209 108 1960 1260];

line_names_femalechase = {
  'GMR_26F09_AE_01'
  'GMR_26E01_AE_01'
  'GMR_21A01_AE_01'
  'GMR_44D11_AE_01'
  'GMR_72C11_AE_01'
  'GMR_20C08_AE_01'
  'GMR_45F11_AE_01'
  'GMR_51B06_AE_01'
  'GMR_71A09_AE_01'
  'GMR_30G01_AE_01'
  'GMR_65G11_AE_01'
  'GMR_48F12_AE_01'
  'GMR_23A07_AE_01'
  'GMR_24D02_AE_01'
  'GMR_45G01_AE_01'
  'GMR_35C10_AE_01'
  'GMR_35C07_AE_01'
  'GMR_31E02_AE_01'
  'GMR_16F02_AE_01'
  };

CachePerLineAnatomyImages(line_name,'imdatafile',imdatafile,'maskfile',maskfile,'cachedir',cachedir);

CacheGroupAnatomyImages(line_names_femalechase,'groupname','female_chase_hits_20131009',...
  'alllinesmatfile',alllinesmatfile,...
  'cachedir',cachedir,...
  'imdatafile',imdatafile,'maskfile',maskfile);

ShowLineAnatomy(line_names_femalechase,'imdatafile',imdatafile,'maskfile',maskfile,...
  'cachedir','/nobackup/branson/AnatomyCacheData20131008',...
  'groupname','female_chase_hits_20131009','figpos',figpos)