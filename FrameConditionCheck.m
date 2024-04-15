function [isallowed,h] = FrameConditionCheck(trx,fly,n,frameconditionparams,haxcurr)

h = zeros(0,2);
if ~exist('haxcurr','var'),
  haxcurr = [];
end
if ~exist('colors','var'),
  colors = [];
end


isallowed = true(1,n);

if isempty(colors)
  colors = lines(7);
end
colori = 1;

for i = 1:2:numel(frameconditionparams)-1,
  
  param = strtrim(frameconditionparams{i});
  val = frameconditionparams{i+1};
  
  plotisstim = contains(param,'stim');
    
  % parse condition type
  m = regexp(param,'^(?<type>(min)|(max)|(equal)|(start)|(end)|(partm?[\d\.]+tom?[\d\.]+[sfp]))_(?<param>.+)$','names','once');
  if isempty(m),
    paramtype = '';
    param1 = param;
  else
    paramtype = m.type;
    param1 = m.param;
  end
  partparams = regexp(paramtype,'^part(?<start>m?[\d\.]+)to(?<end>m?[\d\.]+)(?<units>[sfp])$','names','once');
  if ~isempty(partparams),
    if partparams.start(1) == 'm',
      partparams.start(1) = '-';
    end
    if partparams.end(1) == 'm',
      partparams.end(1) = '-';
    end
    partparams.start = str2double(partparams.start);
    partparams.end = str2double(partparams.end);
    
    % mixed signs are confusing
    if partparams.start < 0,
      assert(partparams.end <= 0);
    end
    if partparams.end > 0,
      assert(partparams.start >= 0);
    end
    paramtype = 'part';
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
      isallowedcurr = [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) >= str2double(val)),false(1,n2)];
    else
      isallowedcurr = [trx(fly).(param1)(1:n1) >= str2double(val),false(1,n2)];
    end
    isallowed = isallowed & isallowedcurr;

  elseif strcmpi(paramtype,'max'),
    if strcmp(param1,'timestamps'),
      isallowedcurr = [((trx(fly).timestamps(1:n1) - trx(fly).timestamps(1)) <= str2double(val)),false(1,n2)];
    else
      isallowedcurr = [trx(fly).(param1)(1:n1) <= str2double(val),false(1,n2)];
    end
    isallowed = isallowed & isallowedcurr;
  elseif strcmpi(paramtype,'equal'),
    isallowedcurr = [trx(fly).(param1)(1:n1) == str2double(val),false(1,n2)];
    isallowed = isallowed & isallowedcurr;
  elseif strcmpi(paramtype,'start'),
    isallowedcurr1 = trx(fly).(param1)(1:n1) == str2double(val);
    isallowedcurr = [false,isallowedcurr1(2:n1)&~isallowedcurr1(1:n1-1),false(1,n2)];
    isallowed = isallowed & isallowedcurr;
  elseif strcmpi(paramtype,'end'),
    isallowedcurr1 = trx(fly).(param1)(1:n1) == str2double(val);
    isallowedcurr = [isallowedcurr1(1:n1-1)&~isallowedcurr1(2:n1),false(1,n2+1)];
    isallowed = isallowed & isallowedcurr;
  elseif strcmpi(paramtype,'part'),
    isallowedcurr = trx(fly).(param1)(1:n1) == str2double(val);
    [t0s,t1s] = get_interval_ends(isallowedcurr);
    t1s = t1s-1;
    
    if ~isempty(haxcurr) && plotisstim && ~isempty(t0s),
      n3 = numel(t0s);
      ylim = get(haxcurr,'YLim');
      h(colori,1) = patch(trx(fly).firstframe-1+[t0s(:),t1s(:),t1s(:),t0s(:),t0s(:)]',ylim(1+[zeros(n3,2),ones(n3,2),zeros(n3,1)])',colors(colori,:)*.3+.7,'LineStyle','none','Parent',haxcurr);
    end

    if partparams.end == 0,
      signend = sign(partparams.start);
    else
      signend = sign(partparams.end);
    end
    if partparams.start == 0,
      signstart = sign(partparams.end);
    else
      signstart = sign(partparams.start);
    end

    switch partparams.units,
        
      case 'p', % percent
        assert(signend>=0 && signstart >=0);
        for j = 1:numel(t0s),
          ntmp = t1s(j)-t0s(j)+1;
          t0 = t0s(j)+round((ntmp-1)*partparams.start/100);
          t1 = t0s(j)+round((ntmp-1)*partparams.end/100);
          isallowedcurr([t0s(j):t0-1,t1+1:t1s(j)]) = false;
        end
      case 's', % seconds
        for j = 1:numel(t0s),
          ts = trx(fly).timestamps(t0s(j):t1s(j));
          if signstart < 0,
            t0 = ts(end);
          else
            t0 = ts(1);
          end
          ok1 = ts - t0 >= partparams.start;
          if signend < 0,
            t1 = ts(end);
          else
            t1 = ts(1);
          end
          ok2 = ts - t1 <= partparams.end;
          isallowedcurr(t0s(j):t1s(j)) = ok1 & ok2;
%           dts = trx(fly).timestamps(t0s(j):t1s(j))-trx(fly).timestamps(t0s(j));
%           isallowedcurr(t0s(j):t1s(j)) = dts >= partparams.start & dts <= partparams.end;
        end
      case 'f', % frames
        for j = 1:numel(t0s),
          if signstart < 0,
            t0 = max(t0s(j),t1s(j)+partparams.start);
          else
            t0 = min(t1s(j),t0s(j)+partparams.start);
          end
          if signend < 0,
            t1 = max(t0s(j),t1s(j)+partparams.end);
          else
            t1 = min(t1s(j),t0s(j)+partparams.end);
          end
          isallowedcurr([t0s(j):t0-1,t1+1:t1s(j)]) = false;
        end   
    end
    isallowed = isallowed & [isallowedcurr(1:n1),false(1,n2)];
  else
    valn = str2double(val);
    if isnan(valn),
      isallowedcurr = [strcmp(trx(fly).(param1)(1:n1),val),false(1,n2)];
    else
      isallowedcurr = [trx(fly).(param1)(1:n1) == valn,false(1,n2)];
    end
    isallowed = isallowed & isallowedcurr;
  end
  if ~isempty(haxcurr) && plotisstim && any(isallowedcurr)
    [t0s,t1s] = get_interval_ends(isallowedcurr); t1s = t1s-1;
    n3 = numel(t0s);
    ylim = get(haxcurr,'YLim');
    h(colori,2) = patch(trx(fly).firstframe-1+[t0s(:),t1s(:),t1s(:),t0s(:),t0s(:)]',ylim(1+[zeros(n3,2),ones(n3,2),zeros(n3,1)])',colors(colori,:)*.5+.5,'LineStyle','none','Parent',haxcurr);
    colori = colori + 1;
    if colori > size(colors,1), colori = 1; end
  end

  
end

h = h(:)';