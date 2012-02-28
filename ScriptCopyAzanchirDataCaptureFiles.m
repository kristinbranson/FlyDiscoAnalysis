%% soft-link to data capture files, create files containing lists of experiment directories

rootoutputdirs = {'/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_CS_20120204',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_nicotine_mathias_berlin_20120211'};

expdirnames = {'expdirs_housing_pBDPGAL4U_20111216.txt',...
  'expdirs_housing_CS_20120204.txt',...
  'expdirs_mating_galit_CS_20120211.txt',...
  'expdirs_nicotine_mathias_berlin_20120211.txt'};

for i = 1:numel(rootoutputdirs),
  expdirs = importdata(fullfile(rootoutputdirs{i},expdirnames{i}));
  for j = 1:numel(expdirs),
    SymbolicCopyFlyBowlDataCaptureFiles(expdirs{j},fullfile(rootoutputdirs{i},'results'));
  end
  fid = fopen(expdirnames{i},'w');
  for j = 1:numel(expdirs),
    [~,basename] = fileparts(expdirs{j});
    fprintf(fid,'%s\n',fullfile(rootoutputdirs{i},'results',basename));
  end
  fclose(fid);
end

%% soft-link to experiments already processed in 20120109_non_olympiad_azanchir

rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216/results';
rootinputdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120109_non_olympiad_azanchir/ctrax_test';
tmp = dir(fullfile(rootinputdir,'p*_2011*'));
for i = 1:numel(tmp),
  cmd = sprintf('ln -s %s %s',fullfile(rootinputdir,tmp(i).name),fullfile(rootoutputdir,tmp(i).name));
  unix(cmd);
end

%% fix metadata files for housing_CS

tmpdata = SAGEListBowlExperiments('daterange',{'20120204T000000','20120206T000000'});
idx = strcmpi({tmpdata.protocol},'EP_flybowl_v011p2.xls');

% copy over metadata files so we can edit them
rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_CS_20120204';
for i = find(idx),
  inexpdir = tmpdata(i).file_system_path;
  [~,experiment_name] = fileparts(inexpdir);
  expdir = fullfile(rootoutputdir,'results',experiment_name);
  metadatafile = fullfile(expdir,'Metadata.xml');
  inmetadatafile = fullfile(inexpdir,'Metadata.xml');
  % delete the soft linked metadata file
  cmd = sprintf('rm %s',metadatafile);
  unix(cmd);
  % copy the file over
  cmd = sprintf('cp %s %s',inmetadatafile,metadatafile);
  unix(cmd);
  % replace gender
  cmd = ['sed -i ''s/gender=\"m\"/gender=\"b\"/g'' ',metadatafile];
  unix(cmd);
end