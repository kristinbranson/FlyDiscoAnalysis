rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
tmp = dir(fullfile(rootdatadir,'*20*'));
expdirs = cellfun(@(x) fullfile(rootdatadir,x),{tmp.name},'UniformOutput',false);
hasmetadata = false(1,numel(expdirs));
isprimary = true(1,numel(expdirs));
isok = false(1,numel(expdirs));
msgs = cell(1,numel(expdirs));
isscreentype = true(1,numel(expdirs));
for i = i:numel(expdirs),
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
  elseif ~strcmp(metadata.screen_type,'primary'),
    isprimary(i) = false;
    fprintf('%s: not primary\n',name);
    continue;
  end
  [isok(i),msgs{i}] = FlyBowlAutomaticChecks_Complete(expdir,'analysis_protocol','tmp','debug',true);
  fprintf('(%d) %s: OK\n',i,name);
end

fid = fopen('PrimaryDataCheck20130306Msgs.txt','w');

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