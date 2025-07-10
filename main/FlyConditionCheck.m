function isallowed = FlyConditionCheck(trx,fly,flyconditionparams)

isallowed = true;

for i = 1:2:numel(flyconditionparams)-1,
  
  param = strtrim(flyconditionparams{i});
  val = flyconditionparams{i+1};
  
  % parse condition type
  if ~isempty(regexp(param,'^min_','once')),
    isallowed = trx(fly).(param(5:end)) >= str2double(val);
  elseif ~isempty(regexp(param,'^max_','once')),
    isallowed = trx(fly).(param(5:end)) <= str2double(val);
  else
    valn = str2double(val);
    if isnan(valn),
      isallowed = strcmp(trx(fly).(param),val);
    else
      isallowed = trx(fly).(param) == valn;
    end
  end
  if ~isallowed,
    break;
  end
  
end