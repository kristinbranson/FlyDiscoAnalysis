function examinestats = ReadExamineExperimentVariables(examinefile)

if ~exist(examinefile,'file'),
  error('ExamineExperimentVariables file %s does not exist',examinefile);
end
fid = fopen(examinefile,'r');

examinestats = {};
while true,
  s = fgetl(fid);
  if ~ischar(s), break; end
  examinestats{end+1} = regexprep(s,',','_'); %#ok<AGROW>
  %fns = strtrim(fns);
  %examinestats{end+1} = fns; %#ok<AGROW>
end

fclose(fid);