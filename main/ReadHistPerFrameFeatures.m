function hist_params = ReadHistPerFrameFeatures(histparamsfile)

if ~exist(histparamsfile,'file'),
  error('Hist parameter file %s does not exist',histparamsfile);
end

fid = fopen(histparamsfile,'r');
if fid < 0,
  error('Could not open parameter file %s for reading',histparamsfile);
end

fields = {};
frameconditions = {};
flyconditions = {};
minZ = {};
nhist = 0;

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
    error('Hist params should be fieldname,frameconditionname,flyconditionname,minZ. Read:\n%s',s);
  end
  
  nhist = nhist + 1;
  
  fields{nhist} = m{1}; %#ok<AGROW>
  frameconditions{nhist} = m{2}; %#ok<AGROW>
  flyconditions{nhist} = m{3}; %#ok<AGROW>
  minZ{nhist} = str2double(m{4}); %#ok<AGROW>
  
end

hist_params = struct('field',fields,'framecondition',frameconditions,'flycondition',flyconditions,'minZ',minZ);