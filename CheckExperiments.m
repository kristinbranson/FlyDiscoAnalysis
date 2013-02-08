function CheckExperiments(expdirs,varargin)

[checkprotocol,checkfix,classifierparamsfiles,npar,analysis_protocol,usesage] = myparse(varargin,...
  'checkprotocol',true,'checkfix',true,...
  'classifierparamsfiles',{},...
  'npar',[],...
  'analysis_protocol','current',...
  'usesage',true);

if isempty(npar),
  npar = max(1,matlabpool('size'));
end
nexps = numel(expdirs);

experiment_names = cell(1,nexps);
for i = 1:nexps,
  [~,name] = fileparts(expdirs{i});
  experiment_names{i} = ['FlyBowl_',name];
end
if usesage,
  metadata = SAGEListBowlExperiments('experiment_name',experiment_names,'checkflags',false,...
    'data_type',{'sexclassifier_diagnostics_mean_nflies','QuickStats_BackSubStats_meanNConnComps'},...
    'removemissingdata',false);
else
  for i = 1:nexps,
    expdir = expdirs{i};
    metadatacurr = parseExpDir(expdir,true);
    metadatacurr.automated_pf = 'U';
    metadatacurr.manual_pf = 'U';
    if ~exist(expdir,'dir'),
      fprintf('Experiment directory %s does not exist\n',expdir);
      metadatacurr.automated_pf = 'F';
    else
      filename = fullfile(expdir,'automatic_checks_complete_results.txt');
      if ~exist(filename,'file'),
        filename = fullfile(expdir,'automatic_checks_incoming_results.txt');
      end
      if exist(filename,'file'),
        res = ReadParams(filename);
        if isfield(res,'automated_pf'),
          metadatacurr.automated_pf = res.automated_pf;
        end
        if isfield(res,'automated_pf_category'),
          metadatacurr.automated_pf_category = res.automated_pf_category;
        end
      end
    end
    metadata(i) = metadatacurr;
  end
end

[ism,idx] = ismember(experiment_names,{metadata.experiment_name});
if any(~ism),
  fprintf('The following experiments were not found in SAGE, and will be ignored in the rest of the checks:\n');
  fprintf('%s\n',expdirs{~ism});
  expdirs(~ism) = [];
  idx(~ism) = [];  
  nexps = numel(expdirs);
end
[~,order] = sort(idx);
metadata = metadata(order);
for i = 1:numel(metadata),
  metadata(i).file_system_path = expdirs{i};
end

if checkprotocol,
  maxareas = nan(1,nexps);
  analysis_protocols = cell(1,nexps);
  ctrax_versions = cell(1,nexps);
  for i = 1:nexps,
    fn = fullfile(expdirs{i},'analysis_protocol.txt');
    if ~exist(fn,'file'),
      analysis_protocols{i} = 'MISSING';
    else
      fid = fopen(fn,'r');
      s = fgetl(fid);
      s = strtrim(s);
      fclose(fid);
      analysis_protocols{i} = s;
    end
    annfile = fullfile(expdirs{i},'movie.ufmf.ann');
    if ~exist(annfile,'file'),
      ctrax_versions{i} = 'MISSING';
    else
      [maxarea,version] = read_ann(annfile,'maxarea','version');
      maxareas(i) = maxarea;
      ctrax_versions{i} = version;
    end
  end
  
  % check analysis protocols
  [allprotocols,~,protocolidx] = unique(analysis_protocols);
  [allmaxareas,~,maxareaidx] = unique(maxareas);
  [allctrax_versions,~,ctrax_versionidx] = unique(ctrax_versions);
  
  fprintf('Analysis protocols:\n');
  ispass = strcmp({metadata.automated_pf},'P') & ~strcmp({metadata.manual_pf},'F');
  isunknown = strcmp({metadata.automated_pf},'U') & ~strcmp({metadata.manual_pf},'F');
  isfailure = ~(ispass|isunknown);
  
  for i = 1:numel(allprotocols),
    fprintf('%s: %d experiments, %d pass, %d unknown, %d fail\n',allprotocols{i},nnz(protocolidx==i),...
      nnz(ispass&protocolidx==i),nnz(isunknown&protocolidx==i),nnz(isfailure&protocolidx==i));
  end
    
  fprintf('Max areas:\n');
  for i = 1:numel(allmaxareas),
    fprintf('%f: %d experiments, %d pass, %d unknown, %d fail\n',allmaxareas(i),nnz(maxareaidx==i),...
      nnz(ispass&maxareaidx==i),nnz(isunknown&maxareaidx==i),nnz(isfailure&maxareaidx==i));
  end

  fprintf('Ctrax versions:\n');
  for i = 1:numel(allctrax_versions),
    fprintf('%s: %d experiments, %d pass, %d unknown, %d fail\n',allctrax_versions{i},nnz(ctrax_versionidx==i),...
      nnz(ispass&ctrax_versionidx==i),nnz(isunknown&ctrax_versionidx==i),nnz(isfailure&ctrax_versionidx==i));
  end
  
  input('Hit enter to continue.');
  
end

if checkfix,
  
  issuccess = nan(1,nexps);
  msgs = cell(1,numel(expdirs));
  nproblems = 0;
  ii = 0;
  datenums = datenum({metadata.exp_datetime},'yyyymmddTHHMMSS');

  hfig = 1;
  figure(hfig);
  clf;

  hdata = plot(datenums,issuccess,'k.');
  set(gca,'XLim',[min(datenums)-1,max(datenums)+1],'YLim',[-.25,1.25]);
  datetick('x','mmmyy','keeplimits');
  hti = title(sprintf('%d problems / %d checked',nproblems,ii));
  order = 1:nexps;

  for ii = 1:npar:nexps

    
    %i = order(ii);
  
    idxcurr = order(ii:min(ii+npar-1,numel(metadata)));
    currissuccess = nan(1,numel(idxcurr));
    currmsgs = cell(1,numel(idxcurr));
    curr_expdirs = {metadata(idxcurr).file_system_path};
    
    parfor j = 1:numel(idxcurr),
      fprintf('***** Checking experiment %s...\n',curr_expdirs{j});
      [currissuccess(j),currmsgs{j}] = CheckExperiment20120726(curr_expdirs{j},'analysis_protocol',analysis_protocol);
    end
  
    issuccess(idxcurr) = currissuccess;
    msgs(idxcurr) = currmsgs;
    nproblems = nnz(issuccess==0);
    
    set(hdata,'YData',issuccess);
    set(hti,'String',sprintf('%d problems / %d checked',nproblems,ii+numel(idxcurr)-1));
  
    drawnow;
  end

  for i = find(~issuccess),
    fprintf('\n%s\n',metadata(i).file_system_path);
    fprintf('%s\n',msgs{i}{:});
  end
  
end

if ~isempty(classifierparamsfiles),
  [issuccess,msgs] = CheckScores(expdirs,classifierparamsfiles);
  issuccess = all(issuccess,1);
  if any(~issuccess),
    fprintf('Problems with scores:\n');
    for i = find(~issuccess),
      fprintf('\n%s:\n',expdirs{i});
      fprintf('%s\n',msgs{i}{:});
    end
  end
end