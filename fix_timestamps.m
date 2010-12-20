function fix_timestamps(filename)

load(filename,'trx');
nflies = length(trx);  %#ok<NODEF>
nframes = max([trx.endframe]);
timestamps = nan(1,nframes);
for fly = 1:nflies,
  if numel(trx(fly).timestamps) ~= trx(fly).nframes,
    continue;
  end
  idx = trx(fly).firstframe:trx(fly).endframe;
  if any(~isnan(timestamps(idx)) & timestamps(idx) ~= trx(fly).timestamps),
    error('Mismatch for fly %d',fly);
  end
  timestamps(idx) = trx(fly).timestamps;
end

isok = true;
for fly = 1:nflies,
  if numel(trx(fly).timestamps) ~= trx(fly).nframes || ...
      ~all(trx(fly).timestamps == timestamps(trx(fly).firstframe:trx(fly).endframe)),
    isok = false;
    break;
  end
end

if isok,
  fprintf('No problems with %s\n',filename);
else
  fprintf('Fixing file %s...\n',filename);
end

if any(isnan(timestamps)),
  error('Did not set all timestamps');
end

if any(timestamps(1:end-1) >= timestamps(2:end)),
  error('Timestamps not always increasing');
end

for fly = 1:nflies,
  trx(fly).timestamps = timestamps(trx(fly).firstframe:trx(fly).endframe); %#ok<AGROW>
end

for fly = 1:nflies,
  if numel(trx(fly).timestamps) ~= trx(fly).nframes,
    error('Timestamps still mismatched for fly %d',fly);
  end
end

copyfile(filename,[filename,'.bak']);
save('-append',filename,'trx','timestamps');