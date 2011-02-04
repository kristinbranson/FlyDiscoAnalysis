function fix_timestamps(trxfile,ctraxfile)

load(trxfile,'trx');
load(ctraxfile,'timestamps');

trx = apply_convert_units(trx); %#ok<NODEF>

nflies = length(trx);

isok = true;
for fly = 1:nflies,
  if numel(trx(fly).timestamps) ~= trx(fly).nframes || ...
      ~all(trx(fly).timestamps == timestamps(trx(fly).firstframe:trx(fly).endframe)),
    isok = false;
    break;
  end
end

if isok,
  fprintf('No problems with %s\n',trxfile);
else
  fprintf('Fixing file %s...\n',trxfile);
end

if any(isnan(timestamps)),
  error('Did not set all timestamps');
end

if any(timestamps(1:end-1) >= timestamps(2:end)), %#ok<COLND>
  error('Timestamps not always increasing');
end

for fly = 1:nflies,
  trx(fly).timestamps = timestamps(trx(fly).firstframe:trx(fly).endframe);
end

for fly = 1:nflies,
  if numel(trx(fly).timestamps) ~= trx(fly).nframes,
    error('Timestamps still mismatched for fly %d',fly);
  end
end

copyfile(trxfile,[trxfile,'.bak']);
save('-append',trxfile,'trx','timestamps');