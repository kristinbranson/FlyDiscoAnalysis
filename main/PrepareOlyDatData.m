function data = PrepareOlyDatData(data)

if isempty(data),
  return;
end

fns = fieldnames(data);
for i = 1:numel(fns),
  fn = fns{i};
  if isstruct(data(1).(fn)),
    for j = 1:numel(data),
      [newfns,vals] = FlattenStructs(data(j).(fn),fn);
      for k = 1:numel(newfns),
        data(j).(newfns{k}) = vals{k};
      end
    end
    data = rmfield(data,fn);
  end
end

% add extra statistics
for i = 1:numel(data),
  
  % numeric date
  data(i).edt = datenum(data(i).exp_datetime,'yyyymmddTHHMMSS');
  
  % rig x bowl
  if isfield(data,'rig') && isfield(data,'bowl'),
    data(i).rig_bowl = sprintf('%d%s',data(i).rig,data(i).bowl);
  end
  
  % numerical date, rounded to day
  data(i).edt_day = floor(data(i).edt);
  
  % numerical date, rounded to week
  data(i).edt_week = floor((data(i).edt-2)/7)*7+2;
  
  % numerical date, rounded to month
  data(i).edt_month = datenum(datestr(data(i).edt,'yyyymm'),'yyyymm');
  
  % part of day
  hrofday = data(i).edt - floor(data(i).edt);
  data(i).part_of_day = double(hrofday>=.5);
  
  % hour of day
  data(i).hour_of_day = floor(hrofday*24);
  
end

function [fns,vals] = FlattenStructs(in,prefix)

fns = {};
vals = {};

fns0 = fieldnames(in);
for i = 1:numel(fns0),
  fn = fns0{i};
  if ~isempty(prefix),
    fn1 = [prefix,'_',fn];
  end
  if ~isstruct(in.(fn)),
    fns{end+1} = fn1; %#ok<AGROW>
    vals{end+1} = in.(fn); %#ok<AGROW>
  else
    [newfns,newvals] = FlattenStructs(in.(fn),fn1);
    newfns = strrep(newfns,'flyany_frame','');
    newfns = strrep(newfns,'stats_perframe','stats');
    fns(end+1:end+numel(newfns)) = newfns;
    %fns(end+1:end+numel(newfns)) = cellfun(@(s) [fn,'_',s],newfns,'UniformOutput',false);
    vals(end+1:end+numel(newvals)) = newvals;
  end
end