function examinestats = ReadExamineBehaviorVariables(examinefile)

if ~exist(examinefile,'file'),
  error('ExamineExperimentVariables file %s does not exist',examinefile);
end
fid = fopen(examinefile,'r');

examinestats = {};
while true,
  s = fgetl(fid);
  if ~ischar(s), break; end
  if isempty(s),
    continue;
  end
  if ~isempty(regexp(s,'^\s*#','once')),
    continue;
  end
  fns = regexp(s,',','split');
  %examinestats{end+1} = regexprep(s,',','_'); %#ok<AGROW>
  fns = strtrim(fns);
  examinestats{end+1} = fns; %#ok<AGROW>
end

fclose(fid);