function AddSetToData(datafile)

MAX_SET_TIMERANGE = 10/(24*60);
MAX_EXPS_PER_SET = 4;

load(datafile,'rawdata');

% add "set" -- super-experiment

[line_names,~,lineidx] = unique({rawdata.line_name});
sets = nan(1,numel(rawdata));
seti = 0;
for linei = 1:numel(line_names),
  expidx1 = find(lineidx==linei);
  [rigs,~,rigidx] = unique([rawdata(expidx1).rig]);
  for rigi = 1:numel(rigs),
    expidx2 = expidx1(rigidx==rigi);

    % sort by datetime
    [exp_datenum,order] = sort(datenum({rawdata(expidx2).exp_datetime},'yyyymmddTHHMMSS'));
    expidx2 = expidx2(order);
    min_set_time = 0;
    nperset = 0;
  
    for i = 1:numel(expidx2),
      
      % start a new set?
      if exp_datenum(i)-min_set_time > MAX_SET_TIMERANGE || nperset >= MAX_EXPS_PER_SET,
        min_set_time = inf;
        seti = seti+1;
        nperset = 0;
      end
      sets(expidx2(i)) = seti;
      nperset = nperset+1;
      min_set_time = min(min_set_time,exp_datenum(i));
      
    end
  end
end

for seti = 1:max(sets),
  expidx = find(sets==seti);
  %fprintf('Experiments in set %d:\n',seti);
  %fprintf('%s\n',rawdata(expidx).experiment_name);
  [~,order] = sort({rawdata(expidx).exp_datetime});
  min_datetime = rawdata(expidx(order(1))).exp_datetime;
  set_name = sprintf('%s__Rig%d__%s',rawdata(expidx(1)).line_name,...
    rawdata(expidx(1)).rig,...
    min_datetime);
  for i = expidx(:)',
    rawdata(i).set = set_name;
  end
end

save('-append',datafile,'rawdata');