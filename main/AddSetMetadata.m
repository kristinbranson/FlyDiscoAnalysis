function metadata = AddSetMetadata(metadata,varargin)

[MAX_SET_TIMERANGE,MAX_EXPS_PER_SET] = myparse(varargin,...
  'MAX_SET_TIMERANGE',10/(24*60),...
  'MAX_EXPS_PER_SET',4);

if ~isfield(metadata,'genotype'),
  for i = 1:numel(metadata),
    metadata(i).genotype = sprintf('%s__%s',metadata(i).line_name,metadata(i).effector);
  end
end

[genotypes,~,genotypeidx] = unique({metadata.genotype});
sets = nan(1,numel(metadata));
seti = 0;
for geni = 1:numel(genotypes),
  expidx1 = find(genotypeidx==geni);
  [rigs,~,rigidx] = unique([metadata(expidx1).rig]);
  for rigi = 1:numel(rigs),
    expidx2 = expidx1(rigidx==rigi);
    
    % sort by datetime
    [exp_datenum,order] = sort(datenum({metadata(expidx2).exp_datetime},'yyyymmddTHHMMSS'));
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
  %fprintf('%s\n',metadata(expidx).experiment_name);
  [~,order] = sort({metadata(expidx).exp_datetime});
  min_datetime = metadata(expidx(order(1))).exp_datetime;
  set_name = sprintf('%s__Rig%d__%s',metadata(expidx(1)).genotype,...
    metadata(expidx(1)).rig,...
    min_datetime);
  for i = expidx(:)',
    metadata(i).set = set_name;
  end
end