data = SAGEListBowlExperiments('removemissingdata',false,'checkflags',false);
for i = 1:numel(data),
  if mod(i,100) == 0,
    fprintf('%d / %d\n',i,numel(data));
  end
  filename = fullfile(data(i).file_system_path,'analysis_protocol.txt');
  didread = false;
  if exist(filename,'file'),
    try
      fid = fopen(filename,'r');
      data(i).full_analysis_protocol = fgetl(fid);
      didread = true;
      fclose(fid);
    catch ME
      error(getReport(ME));
    end
  end
  if ~didread,
    data(i).full_analysis_protocol = data(i).analysis_protocol;
  end
end

data0 = data;

%% get stats about each protocol

fns = {'ctrax_version','ctrax_settings','analysis_settings','analysis_version'};

% remove failures
goodidx = ~strcmpi('F',{data0.automated_pf}) & ~strcmpi('F',{data0.manual_pf});
data = data0(goodidx);

% all analysis protocols
analysis_protocols = unique({data.full_analysis_protocol});

% initialize
protocolinfo = struct;
for i = 1:numel(analysis_protocols),  
  
  s = regexp(analysis_protocols{i},',','split');
  ss = regexp(s,':','split');
  names = cellfun(@(x)x{1},ss,'UniformOutput',false);
  [ism,idx] = ismember(fns,names);
  issuccess = true;
  for k = 1:numel(fns),
    if ism(k),
      protocolinfo(i).(fns{k}) = ss{idx(k)}{2};
    else
      protocolinfo(i).(fns{k}) = 'unknown';
      issuccess = false;
    end
  end
  
  if issuccess,
    protocolinfo(i).name = sprintf('%s_%s_%s',protocolinfo(i).ctrax_version,...
      protocolinfo(i).ctrax_settings,protocolinfo(i).analysis_settings);
  else
    protocolinfo(i).name = analysis_protocols{i};
  end
  protocolinfo(i).fullnames = analysis_protocols{i};
  protocolinfo(i).dates = [];
  protocolinfo(i).screen_types = {};
  protocolinfo(i).nexps = 0;
  protocolinfo(i).nexps_primary = 0;
  protocolinfo(i).mindate = 'NA';
  protocolinfo(i).maxdate = 'NA';
  protocolinfo(i).mindate_primary = 'NA';
  protocolinfo(i).maxdate_primary = 'NA';
end

[analysis_protocols,~,idx] = unique({protocolinfo.name});
dokeep = false(1,numel(analysis_protocols));
for i = 1:numel(analysis_protocols),
  idx1 = find(idx==i);
  keepidx = idx1(1);
  removeidx = idx1(2:end);
  if ~isempty(removeidx),
    protocolinfo(keepidx).analysis_version = {protocolinfo(keepidx).analysis_version,...
      protocolinfo(removeidx).analysis_version};
    protocolinfo(keepidx).fullnames = {protocolinfo(keepidx).fullnames,...
      protocolinfo(removeidx).fullnames};
  else
    protocolinfo(keepidx).analysis_version = {protocolinfo(keepidx).analysis_version};
    protocolinfo(keepidx).fullnames= {protocolinfo(keepidx).fullnames};
  end
  dokeep(keepidx) = true;    
end
protocolinfo = protocolinfo(dokeep);
analysis_protocols = {protocolinfo.name};


% loop through the unique protocols
for i = 1:numel(analysis_protocols),
  
  
  % experiments for this protocol
  idx = ismember({data.full_analysis_protocol},protocolinfo(i).fullnames);
  if ~any(idx),
    continue;
  end
  
  % dates for this protocol
  protocolinfo(i).dates = datenum({data(idx).exp_datetime},'yyyymmddTHHMMSS'); 
  protocolinfo(i).mindate = datestr(min(protocolinfo(i).dates),'yyyymmdd');
  protocolinfo(i).maxdate = datestr(max(protocolinfo(i).dates),'yyyymmdd');
  
  % screen types for this protocol
  protocolinfo(i).screen_types = unique({data(idx).screen_type});
  protocolinfo(i).nexps = nnz(idx);

  % primary data
  idx_primary = strcmpi({data(idx).screen_type},'primary');
  protocolinfo(i).nexps_primary = nnz(idx_primary);
  if protocolinfo(i).nexps_primary > 0,
    protocolinfo(i).mindate_primary = datestr(min(protocolinfo(i).dates(idx_primary)),'yyyymmdd');
    protocolinfo(i).maxdate_primary = datestr(max(protocolinfo(i).dates(idx_primary)),'yyyymmdd');
  else
    protocolinfo(i).mindate_primary = 'NA';
    protocolinfo(i).maxdate_primary = 'NA';
  end
  fprintf('%s: %s to %s\n',protocolinfo(i).name,protocolinfo(i).mindate,protocolinfo(i).maxdate);
end

%% histogram dates per protocol


datecenters = unique(floor(datenum({data.exp_datetime},'yyyymmddTHHMMSS')))+.5;
nbins = numel(datecenters);
%nbins = 100;
%mindate = min(cat(1,protocolinfo.dates));
%maxdate = max(cat(1,protocolinfo.dates));
%dateedges = linspace(mindate,maxdate,nbins+1);
%datecenters = (dateedges(1:end-1)+dateedges(2:end))/2;

counts = nan(numel(analysis_protocols),nbins);
for i = 1:numel(analysis_protocols),
  if protocolinfo(i).nexps == 0,
    continue;
  end
  counts(i,:) = hist(protocolinfo(i).dates,datecenters);
end


%counts = bsxfun(@rdivide,counts,sum(counts,2));
h = plot(datecenters,counts,'.-');
colors = jet(numel(h))*.7;
markers = {'.','x','o','*','s','d'};
for i = 1:numel(h),
  set(h(i),'Color',colors(i,:),'Marker',markers{mod(i-1,numel(markers))+1});
end
hleg = legend(analysis_protocols);
set(hleg,'interpreter','none');
xtick = datecenters(1)-.5:14:datecenters(end)-.5;
set(gca,'XTick',xtick);
axisalmosttight;
datetick('x','yymmdd','keepticks','keeplimits');

%% print out timeline

[~,protocoli] = max(counts,[],1);
ischange = [true,protocoli(1:end-1)~=protocoli(2:end)];
idx = [find(ischange),nbins+1];
for ii = 1:numel(idx)-1,
  i1 = idx(ii);
  i2 = idx(ii+1)-1;
  fprintf('%s: %s-%s (%d days)\n',analysis_protocols{protocoli(i1)},datestr(datecenters(i1),'yyyymmdd'),datestr(datecenters(i2),'yyyymmdd'),datecenters(i2)-datecenters(i1)+1);
end


%% print out table


[~,order] = sort({protocolinfo.maxdate});

for i = order(:)',
  fprintf('| %s\t| %s\t| %s\t| ',protocolinfo(i).ctrax_version,protocolinfo(i).ctrax_settings,protocolinfo(i).analysis_settings);
  fprintf('%s ',protocolinfo(i).analysis_version{:});
  fprintf('\t| %s\t| %s\t| %s\t| %s\t| %d\t| %d\t|',...
    protocolinfo(i).mindate,protocolinfo(i).maxdate,...
    protocolinfo(i).mindate_primary,protocolinfo(i).maxdate_primary,...
    protocolinfo(i).nexps,protocolinfo(i).nexps_primary);
  fprintf(' %s',protocolinfo(i).screen_types{:});
  fprintf('\t| %s\t|\n',analysis_protocols{i});
end

%% compare

protocol1 = '20110407';
protocol2 = '20110426';
ann1 = read_ann(fullfile('../FlyBowlCtrax',protocol1,'settings.ann'));
ann2 = read_ann(fullfile('../FlyBowlCtrax',protocol2,'settings.ann'));
fns = union(fieldnames(ann1),fieldnames(ann2));

for i = 1:numel(fns),
  fn = fns{i};
  if ~isfield(ann1,fn),
    fprintf('%s not in %s\n',fn,protocol1);
    continue;
  end
  if ~isfield(ann2,fn),
    fprintf('%s not in %s\n',fn,protocol2);
    continue;
  end
  tmp1 = ann1.(fn);
  tmp2 = ann2.(fn);
  if ~strcmp(class(tmp1),class(tmp2)),
    fprintf('%s: class mismatch\n',fn);
    fprintf('%s: ',protocol1);
    disp(ann1.(fn));
    fprintf('%s: ',protocol2);
    disp(ann2.(fn));

    continue;
  end
  if ndims(tmp1) ~= ndims(tmp2),
    fprintf('%s: ndims mismatch\n',fn);
    fprintf('%s: ',protocol1);
    disp(ann1.(fn));
    fprintf('%s: ',protocol2);
    disp(ann2.(fn));
    continue;
  end
  if ~all(size(tmp1) == size(tmp2)),
    fprintf('%s: size mismatch\n',fn);
    fprintf('%s: ',protocol1);
    disp(ann1.(fn));
    fprintf('%s: ',protocol2);
    disp(ann2.(fn));
    continue;
  end
  if isnumeric(tmp1),
    if ~all(tmp1(:)==tmp2(:)),
      fprintf('%s: value mismatch\n',fn);
      if numel(tmp1) < 10,
        fprintf('%s: ',protocol1);
        disp(ann1.(fn));
        fprintf('%s: ',protocol2);
        disp(ann2.(fn));
      else
        fprintf('%s: ',protocol1);
        disp(ann1.(fn)(1:10));
        fprintf('%s: ',protocol2);
        disp(ann2.(fn)(1:10));
      end
      continue;
    end
  elseif ischar(tmp1),
    if ~strcmp(tmp1,tmp2),
      fprintf('%s: value mismatch\n',fn);
      fprintf('%s: ',protocol1);
      disp(ann1.(fn));
      fprintf('%s: ',protocol2);
      disp(ann2.(fn));
      continue;
    end
  else
    fprintf('didn''t check values for %s\n',fn);
    fprintf('%s: ',protocol1);
    disp(ann1.(fn));
    fprintf('%s: ',protocol2);
    disp(ann2.(fn));
  end
end
    