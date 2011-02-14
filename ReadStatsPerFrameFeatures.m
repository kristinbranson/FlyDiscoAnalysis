function stats_params = ReadStatsPerFrameFeatures(statsparamsfile)

if ~exist(statsparamsfile,'file'),
  error('Stats parameter file %s does not exist',statsparamsfile);
end

fid = fopen(statsparamsfile,'r');
if fid < 0,
  error('Could not open parameter file %s for reading',statsparamsfile);
end

fields = {};
frameconditions = {};
flyconditions = {};
minZ = {};
nstats = 0;

while true,
  s = fgetl(fid);
  
  % end of file
  if ~ischar(s), break; end
  
  % comments
  if isempty(s) || ~isempty(regexp(s,'^\s*$','once')) || ...
        ~isempty(regexp(s,'^\s*#','once')),
    continue;
  end

  % should be: fieldname,frameconditionname,flyconditionname,minZ
  m = regexp(s,',','split');
  if numel(m) ~= 4,
    error('Stats params should be fieldname,frameconditionname,flyconditionname,minZ. Read:\n%s',s);
  end
  
  nstats = nstats + 1;
  
  fields{nstats} = m{1}; %#ok<AGROW>
  frameconditions{nstats} = m{2}; %#ok<AGROW>
  flyconditions{nstats} = m{3}; %#ok<AGROW>
  minZ{nstats} = str2double(m{4}); %#ok<AGROW>
  
end

stats_params = struct('field',fields,'framecondition',frameconditions,'flycondition',flyconditions,'minZ',minZ);