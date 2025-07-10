function quickstats = loadQuickStats(expdir,varargin)

[quickstatsfilestr] = myparse(varargin,'quickstatsfilestr','QuickStats.txt');

quickstatsfilename = fullfile(expdir,quickstatsfilestr);

if ~exist(quickstatsfilename,'file'),
  error('File %s does not exist',quickstatsfilename);
end

fid = fopen(quickstatsfilename,'r');

quickstats = struct;
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if isempty(s) || ~isempty(regexp(s,'\w*#','once')),
    continue;
  end
  
  parts = regexp(s,',','split');
  if isempty(parts),
    warning('Could not parse line %s',s);
    continue;
  end
  fn = parts{1};
  val = str2double(parts(2:end));
  quickstats.(fn) = val;
end

fclose(fid);