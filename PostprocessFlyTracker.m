function [trxout,ninterpframes,sourceid] = PostprocessFlyTracker(trxin,varargin)

maxFlyTrackerNanInterpFrames = 5;

[maxFlyTrackerNanInterpFrames] = myparse(varargin,...
  'maxFlyTrackerNanInterpFrames',maxFlyTrackerNanInterpFrames);

fns_perframe = {'x','y','theta','a','b',...
  'x_mm','y_mm','a_mm','b_mm','theta_mm','timestamps',...
  'xwingl','ywingl','xwingr','ywingr'};
fns_check = {'x','y','theta','a','b'};

trxout = [];
sourceid = [];
ninterpframes = 0;
for i = 1:numel(trxin),
  trxcurr = trxin(i);
  isbad = false(1,trxcurr.nframes);
  for j = 1:numel(fns_check),
    isbad = isbad | isnan(trxcurr.(fns_check{j}));
  end
  trxcurr.id = numel(trxout)+1;
  if ~any(isbad),
    trxout = structappend(trxout,trxcurr);
    sourceid(end+1) = i; %#ok<AGROW>
    continue;
  end
  if all(isbad),
    continue;
  end
  % first frame bad, just crop
  if isbad(1),
    off = find(~isbad,1)-1;
    trxcurr.firstframe = trxcurr.firstframe+off;
    trxcurr.off = 1-trxcurr.firstframe;
    trxcurr.nframes = trxcurr.nframes-off;
    for j = 1:numel(fns_perframe),
      fn = fns_perframe{j};
      trxcurr.(fn) = trxcurr.(fn)(off+1:end);
    end
    isbad = isbad(off+1:end);    
  end
  % last frame bad, just crop
  if isbad(end),
    off = trxcurr.nframes-find(~isbad,1,'last');
    trxcurr.endframe = trxcurr.endframe-off;
    trxcurr.nframes = trxcurr.nframes-off;
    for j = 1:numel(fns_perframe),
      fn = fns_perframe{j};
      trxcurr.(fn) = trxcurr.(fn)(1:end-off);
    end
    isbad = isbad(1:end-off);
  end
  if ~any(isbad),
    trxout = structappend(trxout,trxcurr);
    sourceid(end+1) = i; %#ok<AGROW>
    continue;
  end
  % bad sequences in the middle
  [st,en] = get_interval_ends(isbad);
  badlength = en-st;
  idxinterp = find(badlength <= maxFlyTrackerNanInterpFrames);
  for jj = 1:numel(idxinterp),
    j = idxinterp(jj);
    t0 = st(j);
    t1 = en(j)-1;
    for k = 1:numel(fns_perframe),
      fn = fns_perframe{k};
      x = linspace(trxcurr.(fn)(t0-1),trxcurr.(fn)(t1+1),badlength(j)+2);
      trxcurr.(fn)(t0-1:t1+1) = x;
    end
    isbad(t0:t1) = false;
    ninterpframes = ninterpframes + (t1-t0+1);
  end
  isgood = ~isbad;
  [st,en] = get_interval_ends(isgood);
  for j = 1:numel(st),
    t0 = st(j);
    t1 = en(j)-1;
    trxnew = trxcurr;
    trxnew.firstframe = t0 + trxcurr.firstframe - 1;
    trxnew.endframe = t1 + trxcurr.firstframe - 1;
    trxnew.nframes = t1-t0+1;
    trxnew.off = 1-trxnew.firstframe;
    trxnew.id = numel(trxout)+1;
    for k = 1:numel(fns_perframe),
      fn = fns_perframe{k};
      trxnew.(fn) = trxcurr.(fn)(t0:t1);
    end
    trxout = structappend(trxout,trxnew);
    sourceid(end+1) = i; %#ok<AGROW>
  end
  
end

for i = 1:numel(trxout),
  trxout(i).dt = diff(trxout(i).timestamps); %#ok<AGROW>
end

%% debug

if false,
  
  for i = 1:numel(trxout),
    
    id = sourceid(i);
    off1 = trxout(i).firstframe-trxin(id).firstframe;
    assert(all(isnan(trxin(id).x(off1+1:off1+trxout(i).nframes)) | trxout(i).x == trxin(id).x(off1+1:off1+trxout(i).nframes)));
    assert(numel(trxout(i).x)==trxout(i).nframes);
    assert((trxout(i).endframe-trxout(i).firstframe+1)==trxout(i).nframes);
    
  end
  
  clf;
  hax(1) = subplot(2,1,1);
  hold on;
  for i = 1:numel(trxin),
    plot(trxin(i).firstframe:trxin(i).endframe,trxin(i).y,'.-');
  end
  hax(2) = subplot(2,1,2);
  hold on;
  for i = 1:numel(trxout),
    plot(trxout(i).firstframe:trxout(i).endframe,trxout(i).y,'.-');
  end
  linkaxes(hax);
end
