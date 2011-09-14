function [isconsistent,filesmissing] = CheckDiagnosticsFileSAGEConsistency(outfid,filename,prefix,data,dict)

EPS = 1e-3;

isconsistent = true;
filesmissing = false;

if ~ischar(filename),
  diagnostics = filename;
  success = true;
else
  
  success = false;
  if ~exist(filename,'file'),
    fprintf(outfid,'%s diagnostics file %s does not exist\n',prefix,filename);
    filesmissing = true;
  else
    try
      diagnostics = ReadParams(filename);
      success = true;
    catch ME,
      fprintf(outfid,'Could not read %s diagnostics file %s\n',prefix,filename);
      getReport(ME)
      filesmissing = true;
    end
  end
  
end

if nargin < 5,
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
      fprintf(outfid,'Missing field %s in SAGE\n',sagefn);
      isconsistent = false;
    elseif isempty(data.(sagefn)),
      fprintf(outfid,'Field %s is empty in SAGE\n',sagefn);
      isconsistent = false;
    elseif isnumeric(diagnostics.(fn)),
      nsage = numel(data.(sagefn));
      nfile = numel(diagnostics.(fn));
      if nsage ~= nfile,
        fprintf(outfid,'Sizes of field %s do not match for %s\n',fn,prefix);
        isconsistent = false;
      elseif max(abs(double(data.(sagefn)(:)) - diagnostics.(fn)(:))) > EPS,
        if nfile == 1,
          fprintf(outfid,'Field %s = %f in SAGE, %f in %s\n',sagefn,data.(sagefn),diagnostics.(fn),prefix);
        else
          fprintf(outfid,'Field %s in SAGE does not match file for %s\n',sagefn,prefix);
        end
        isconsistent = false;
      end
    end        
  end
end