function params = ReadStimulusPlotParams(paramsfile)

params = ReadParams(paramsfile);
if ~isfield(params,'stimsets'),
  return;
end
stimsets = cell(numel(params.stimsets),2);
for i = 1:numel(params.stimsets),
  s = params.stimsets{i};
  m = regexp(s,'^([^:]*):([\d_]*)$','tokens','once');
  assert(numel(m)==2);
  stimsets{i,1} = m{1};
  periods = str2double(regexp(m{2},'_','split'));
  assert(~any(isnan(periods)));
  stimsets{i,2} = periods;
end
params.stimsets = stimsets;