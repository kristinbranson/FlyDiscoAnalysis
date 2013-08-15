%%

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';

%%

rootdatadir = '/groups/branson/bransonlab/Galit_data_nolinks2';
classifierparamsfiles = {fullfile(settingsdir,'20130722/JAABA_classifier_params.txt')};
[paths,dirs] = mydir(rootdatadir,'isdir',true);

protocols = {};
expdirs = {};
protocol2expdirs = struct;
for i = 1:numel(dirs),
  protocol = dirs(i).name;
  expdirscurr = mydir(paths{i},'isdir',true,'name','\d{8}T\d{6}');
  if ~isempty(expdirscurr),
    expdirs = [expdirs,expdirscurr];
    protocols{end+1} = protocol;
    protocol2expdirs.(protocol) = expdirscurr;
  end
end

%% check 1

issuccess = CheckExperiments(expdirs,'analysis_protocol','20120220_non_olympiad_azanchir_mating_galit_CS_20120211',...
  'classifierparamsfiles',classifierparamsfiles);

%% check 2


hasmetadata = false(1,numel(expdirs));
isok = false(1,numel(expdirs));
msgs = cell(1,numel(expdirs));
isscreentype = true(1,numel(expdirs));
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  [~,name] = fileparts(expdir);
  metadatafile = fullfile(expdir,'Metadata.xml');
  hasmetadata(i) = exist(metadatafile,'file');
  if ~hasmetadata(i),
    fprintf('%s: no metadata\n',name);
    continue;
  end
  try
    metadata = ReadMetadataFile(metadatafile);
  catch ME,
    warning('Could not read metadata file %s: %s',metadatafile,getReport(ME));
    hasmetadata(i) = false;
    fprintf('%s: bad metadata\n',name);
    continue;
  end
  if ~isfield(metadata,'screen_type'),
    fprintf('%s: no screen_type field\n',name);
    isscreentype(i) = false;
  end
  [isok(i),msgs{i}] = FlyBowlAutomaticChecks_Complete(expdir,'analysis_protocol','tmp2','debug',true);
  
  if ~isok(i),
    fprintf('%s: failed automatic checks complete:\n',name);
    fprintf('%s\n',msgs{i}{:});
    continue;
  end
  
  isok(i) = FlyBowlCheckPerFrameFeatures(expdir,'analysis_protocol','20120220_non_olympiad_azanchir_mating_galit_CS_20120211');
  
  if ~isok(i),
    fprintf('%s: Failed per-frame checks\n',name);
  end
  
  fprintf('(%d) %s: OK\n',i,name);
  
end

fid = fopen('GalitDataCheck20130805Msgs.txt','w');

idx = find(~cellfun(@isempty,msgs));
for i = idx,
  if ~isempty(regexp(expdirs{i},'Rig\dPlate\d+Bowl._notstarted','once')),
    continue;
  end
  if exist(fullfile(expdirs{i},'ABORTED'),'file'),
    continue;
  end
  fprintf(fid,'%s:\n',expdirs{i});
  if exist(fullfile(expdirs{i},'automatic_checks_incoming_results.txt'),'file'),
    incoming = ReadParams(fullfile(expdirs{i},'automatic_checks_incoming_results.txt'));
    if isfield(incoming,'automated_pf') && incoming.automated_pf == 'F',
      fprintf(fid,'  automated_pf_incoming = F\n');
      if isfield(incoming,'automated_pf_category'),
        fprintf(fid,'  automated_pf_incoming_category = %s\n',incoming.automated_pf_category);
      end
      continue;
    end
  end
  if ~exist(fullfile(expdirs{i},'ctrax_results.mat'),'file') && exist(fullfile(expdirs{i},'ctrax_log.txt'),'file'),
    [res1,res2] = unix(sprintf('grep ShortUFMFFileError %s',fullfile(expdirs{i},'ctrax_log.txt')));
    if res1 == 0 && ~isempty(res2),
      fprintf(fid,'  ShortUFMFFileError\n');
      continue;
    end
  end
  fprintf(fid,'  %s\n',msgs{i}{:});
end
fclose(fid);

%%

needsresultsmovie = false(1,numel(expdirs));
for i = find(isok),
  expdir = expdirs{i};
  tmp = dir(fullfile(expdir,'ctrax_results_movie*.avi'));
  needsresultsmovie(i) = isempty(tmp);
end