function isallowed = FrameConditionCheck(trx,fly,n,frameconditionparams)


isallowed = true(1,n);

for i = 1:2:numel(frameconditionparams)-1,
  
  param = strtrim(frameconditionparams{i});
  val = frameconditionparams{i+1};
    
  % parse condition type
  m = regexp(param,'^(?<type>(min)|(max)|(equal)|(start))_(?<param>.+)$','names','once');
  if isempty(m),
    paramtype = '';
    param1 = param;
  else
    paramtype = m.type;
    param1 = m.param;
  end
%   if ~isempty(regexp(param,'^min_','once')) || ~isempty(regexp(param,'^max_','once')) ...
%       || ~isempty(regexp(param,'^equal_','once')),
%     param1 = param(5:end);
%   else
%     param1 = param;
%   end
  n1 = min(n,numel(trx(fly).(param1)));
  n2 = n-n1;
  if strcmpi(paramtype,'min'),
    if strcmp(param1,'timestamps'),
      isallowed = isallowed & ...
        [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) >= str2double(val)),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) >= str2double(val),false(1,n2)];
    end
  elseif strcmpi(paramtype,'max'),
    if strcmp(param1,'timestamps'),
      isallowed = isallowed & ...
        [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) <= str2double(val)),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) <= str2double(val),false(1,n2)];
    end
  elseif strcmpi(paramtype,'equal'),
    isallowed = isallowed & [trx(fly).(param1)(1:n1) == str2double(val),false(1,n2)];
  elseif strcmpi(paramtype,'start'),
    tmp = trx(fly).(param1)(1:n1) == str2double(val);
    isallowed = isallowed & [false,tmp(2:n1)&~tmp(1:n1-1),false(1,n2)];
  elseif strcmpi(paramtype,'end'),
    tmp = trx(fly).(param1)(1:n1) == str2double(val);
    isallowed = isallowed & [tmp(1:n1-1)&~tmp(2:n1),false(1,n2+1)];
  else
    valn = str2double(val);
    if isnan(valn),
      isallowed = isallowed & [strcmp(trx(fly).(param1)(1:n1),val),false(1,n2)];
    else
      isallowed = isallowed & [trx(fly).(param1)(1:n1) == valn,false(1,n2)];
    end
  end
  
end