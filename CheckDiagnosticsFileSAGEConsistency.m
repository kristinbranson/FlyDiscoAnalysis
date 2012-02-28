function [isconsistent,filesmissing,msgs] = CheckDiagnosticsFileSAGEConsistency(filename,prefix,data,dict)

EPS = 1e-3;

isconsistent = true;
filesmissing = false;
msgs = {};

if ~ischar(filename),
  diagnostics = filename;
  success = true;
else
  
  success = false;
  if ~exist(filename,'file'),
    %msgs{end+1} = sprintf('%s diagnostics file %s does not exist',prefix,filename);
    filesmissing = true;
  else
    try
      diagnostics = ReadParams(filename);
      success = true;
    catch ME,
      fprintf('Could not read %s diagnostics file %s\n',prefix,filename);
      getReport(ME)
      filesmissing = true;
    end
  end
  
end

if nargin < 4,
  dict = cell(0,2);
end

if success,
  fns = fieldnames(diagnostics);
  for i = 1:numel(fns),
    fn = fns{i};
    j = find(strcmp(dict(:,1),fn),1);
    if isempty(j),
      sagefn = [prefix,'_',fn];
    else
      sagefn = [prefix,'_',dict{j,2}];
    end
    if ~isfield(data,sagefn),
      msgs{end+1} = sprintf('Missing field %s in SAGE',sagefn); %#ok<*AGROW>
      isconsistent = false;
    elseif isempty(data.(sagefn)),
      msgs{end+1} = sprintf('Field %s is empty in SAGE',sagefn);
      isconsistent = false;
    elseif isnumeric(diagnostics.(fn)),
      nsage = numel(data.(sagefn));
      nfile = numel(diagnostics.(fn));
      if iscell(data.(sagefn)),
        data.(sagefn) = str2double(data.(sagefn));
      end
      if nsage ~= nfile,
        msgs{end+1} = sprintf('Sizes of field %s do not match for %s',fn,prefix);
        isconsistent = false;
      else
        badidx = isinf(diagnostics.(fn)(:)) | isnan(diagnostics.(fn)(:));
        if max(abs(double(data.(sagefn)(~badidx)) - diagnostics.(fn)(~badidx))) > EPS,
          if nfile == 1,
            msgs{end+1} = sprintf('Field %s = %f in SAGE, %f in %s',sagefn,data.(sagefn),diagnostics.(fn),prefix);
          else
            msgs{end+1} = sprintf('Field %s in SAGE does not match file for %s',sagefn,prefix);
          end
          isconsistent = false;
        end
      end
    end        
  end
end


