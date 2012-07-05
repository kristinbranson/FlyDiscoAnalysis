%% soft-link to data capture files, create files containing lists of experiment directories

rootoutputdirs = {'/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_housing_CS_20120204',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_nicotine_mathias_berlin_20120211',...
  '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120228_non_olympiad_azanchir_housing_CS_20120225'};

expdirnames = {'expdirs_housing_pBDPGAL4U_20111216.txt',...
  'expdirs_housing_CS_20120204.txt',...
  'expdirs_mating_galit_CS_20120211.txt',...
  'expdirs_nicotine_mathias_berlin_20120211.txt',...
  'expdirs_housing_CS_20120225.txt'};

analysis_protocols = {'20120220_non_olympiad_azanchir_housing_CS_20120204',...
'20120220_non_olympiad_azanchir_housing_pBDPGAL4U_20111216',...
'20120220_non_olympiad_azanchir_mating_galit_CS_20120211',...
'20120220_non_olympiad_azanchir_nicotine_mathias_berlin_20120211',...
'20120228_non_olympiad_azanchir_housing_CS_20120225'};

STARTI = 5;

% copy over data capture files
for i = STARTI:numel(rootoutputdirs),
  expdirs = importdata(fullfile(rootoutputdirs{i},expdirnames{i}));
  for j = 1:numel(expdirs),
    SymbolicCopyFlyBowlDataCaptureFiles(expdirs{j},fullfile(rootoutputdirs{i},'results'));
  end
end
% output expdirs to file
for i = STARTI:numel(rootoutputdirs),
  expdirs = importdata(fullfile(rootoutputdirs{i},expdirnames{i}));
  fid = fopen(expdirnames{i},'w');
  for j = 1:numel(expdirs),
    [~,basename] = fileparts(expdirs{j});
    fprintf(fid,'%s\n',fullfile(rootoutputdirs{i},'results',basename));
  end
  fclose(fid);
end
% look at load time range
fprintf('Experiment\tLoadtime\tVideoLength\n');
for i = STARTI:numel(rootoutputdirs),
  expdirs = importdata(fullfile(rootoutputdirs{i},expdirnames{i}));
  fliesloaded_times = nan(1,numel(expdirs));
  for j = 1:numel(expdirs),
    metadatafile = fullfile(expdirs{j},'Metadata.xml');
    metadata = ReadMetadataFile(metadatafile);
    fliesloaded_times(j) = metadata.seconds_fliesloaded;
    [~,basename] = fileparts(expdirs{j});
    fprintf('%s\t%d\t%d\n',basename,round(fliesloaded_times(j)),1200-round(fliesloaded_times(j)));
  end
  fprintf('seconds_fliesloaded range: %d to %d\n',round(min(fliesloaded_times)),round(max(fliesloaded_times)));
end

%% what is the range of load times?
for i = STARTI:numel(rootoutputdirs),
  
end

%% command to run on the cluster

fprintf('Run the following on the cluster:\n');
for i = STARTI:numel(rootoutputdirs),
  fprintf('./qsub_FlyBowlAnalysisPipeline.pl %s %s\n',expdirnames{i},analysis_protocols{i});
end