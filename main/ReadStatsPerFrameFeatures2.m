function stats_params = ReadStatsPerFrameFeatures2(statsparamsfile)

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
minZboth = {};
minZfly = {};
norm_frameconditions = {};
norm_flyconditions = {};
norm_minZboth = {};
norm_minZfly = {};
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
  if numel(m) < 5,
    error('Stats params should be fieldname,frameconditionname,flyconditionname,minZboth,minZfly,[norm_frameconditionname,norm_flyconditionname,norm_minZboth,norm_minZfly]. Read:\n%s',s);
  end
  
  nstats = nstats + 1;
  
  fields{nstats} = m{1}; %#ok<AGROW>
  frameconditions{nstats} = m{2}; %#ok<AGROW>
  flyconditions{nstats} = m{3}; %#ok<AGROW>
  minZboth{nstats} = str2double(m{4}); %#ok<AGROW>
  if isnan(minZboth{nstats}),
    error('minZboth is not a number: %s\n',m{4});
  end
  minZfly{nstats} = str2double(m{5}); %#ok<AGROW>
  if isnan(minZfly{nstats}),
    error('minZfly is not a number: %s\n',m{5});
  end
  if minZboth{nstats} > minZfly{nstats},
    error('minZboth must be <= minZfly, Read: %s',s);
  end
  norm_frameconditions{nstats} = ''; %#ok<AGROW>
  norm_flyconditions{nstats} = ''; %#ok<AGROW>
  norm_minZboth{nstats} = ''; %#ok<AGROW>
  norm_minZfly{nstats} = ''; %#ok<AGROW>
  if numel(m) > 5,
    norm_frameconditions{nstats} = m{6}; %#ok<AGROW>
    if numel(m) >= 7,
      norm_flyconditions{nstats} = m{7}; %#ok<AGROW>
    else
      norm_flyconditions{nstats} = flyconditions{nstats}; %#ok<AGROW>
    end
    if numel(m) >= 8,
      norm_minZboth{nstats} = str2double(m{8}); %#ok<AGROW>
      if isnan(norm_minZboth{nstats}),
        error('norm_minZboth is not a number: %s\n',m{8});
      end
    else
      norm_minZboth{nstats} = minZboth{nstats}; %#ok<AGROW>
    end
    if numel(m) >= 9,
      norm_minZfly{nstats} = str2double(m{9}); %#ok<AGROW>
      if isnan(norm_minZfly{nstats}),
        error('norm_minZfly is not a number: %s\n',m{9});
      end
    else
      norm_minZfly{nstats} = minZfly{nstats};  %#ok<AGROW>
    end
  end

  
end

stats_params = struct('field',fields,'framecondition',frameconditions,'flycondition',flyconditions,'minZboth',minZboth,'minZfly',minZfly,...
  'norm_framecondition',norm_frameconditions,'norm_flycondition',norm_flyconditions,'norm_minZboth',norm_minZboth,'norm_minZfly',norm_minZfly);
fclose(fid);