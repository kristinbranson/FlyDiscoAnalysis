function isallowed = FrameConditionCheck(trx,fly,n,frameconditionparams)


isallowed = true(1,n);

for i = 1:2:numel(frameconditionparams)-1,
  
  param = strtrim(frameconditionparams{i});
  val = frameconditionparams{i+1};
    
  % parse condition type
  if ~isempty(regexp(param,'^min_','once')) || ~isempty(regexp(param,'^max_','once')),
    param1 = param(5:end);
  else
    param1 = param;
  end
  n1 = min(n,numel(trx(fly).(param1)));
  n2 = n-n1;
  if ~isempty(regexp(param,'^min_','once'))
    if strcmp(param1,'timestamps'),
      isallowed = isallowed & ...
        [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) >= str2double(val)),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) >= str2double(val),false(1,n2)];
    end
  elseif ~isempty(regexp(param,'^max_','once')),
    if strcmp(param1,'timestamps'),
      isallowed = isallowed & ...
        [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) <= str2double(val)),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) <= str2double(val),false(1,n2)];
    end
  else
    valn = str2double(val);
    if isnan(valn),
      isallowed = isallowed & [strcmp(trx(fly).(param1)(1:n1),val),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) == valn,false(1,n2)];
    end
  end
  
end