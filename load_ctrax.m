function trx = load_ctrax(matName,movieName,annName)

if ~exist(matName,'file'),
  error('Ctrax file %s does not exist',matName);
end

if ~exist(movieName,'var'),
  movieName = '?';
end
if ~exist(annName,'var'),
  annName = '?';
end

in = load(matName);

% always make row vectors
in.x_pos = in.x_pos(:)';
in.y_pos = in.y_pos(:)';
in.maj_ax = in.maj_ax(:)';
in.min_ax = in.min_ax(:)';
in.angle = in.angle(:)';
in.identity = in.identity(:)';
if isfield(in,'timestamp'),
  in.timestamp = in.timestamp(:)';
end

% find identities
idscurr = unique(in.identity);

% matlab uses 1-indexing
in.x_pos = in.x_pos + 1;
in.y_pos = in.y_pos + 1;

% frame number
framenumber = zeros(size(in.x_pos));
j = 0;
for i = 1:length(in.ntargets),
  framenumber(j+(1:in.ntargets(i))) = i;
  j = j + in.ntargets(i);
end;

newidentity = nan(size(in.identity));
for id = idscurr,
  idx = in.identity == id;
  datacurr.x = in.x_pos(idx);
  datacurr.y = in.y_pos(idx);
  datacurr.theta = in.angle(idx);
  datacurr.a = in.maj_ax(idx);
  datacurr.b = in.min_ax(idx);
  if isfield(in,'timestamp'),
    datacurr.timestamp = in.timestamp(idx);
  end
  datacurr.id = id;
  datacurr.moviename = movieName;
  datacurr.annname = annName;
  datacurr.firstframe = framenumber(find(idx,1));
  datacurr.arena.x = nan;
  datacurr.arena.y = nan;
  datacurr.arena.r = nan;
  datacurr.off = -datacurr.firstframe + 1;
  datacurr.nframes = length(datacurr.x);
  datacurr.endframe = datacurr.nframes + datacurr.firstframe - 1;
  if ~exist('trx','var'),
    trx = datacurr;
  else
    trx(end+1) = datacurr; %#ok<AGROW>
  end;
  newidentity(idx) = length(trx);
end
if isfield(trx,'timestamp'),
  for fly = 1:length(trx),
    trx(fly).dt_s = diff(trx(fly).timestamp);
  end
end