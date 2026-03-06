function [data,pos,coltitles,blockids,neuropilnames,neuropilids] = ReadMBVoxelInfo(datafile)

fid = fopen(datafile,'r');
s = fgetl(fid);
coltitles1 = regexp(s,'\t','split');
s = fgetl(fid);
coltitles2 = regexp(s,'\t','split');

coltitles = cell(size(coltitles2));
for i = 1:numel(coltitles2),
  if isempty(coltitles1{i}),
    coltitles{i} = coltitles2{i};
  elseif isempty(coltitles2{i}),
    coltitles{i} = coltitles1{i};
  else
    coltitles{i} = [coltitles1{i},'_',coltitles2{i}];
  end
end

nspecial = 6;
if ~all(strcmpi(coltitles(1:nspecial),{'BlockID','X','Y','Z','Neuropil','Neuropil'})),
  error('Column titles do not match expected values');
end
goodidx = ~cellfun(@isempty,coltitles);
goodidx(1:nspecial) = false;

coltitles = coltitles(nspecial+1:end);
blockids = nan(0,1);
pos = nan(0,3);
neuropilnames = cell(0,1);
neuropilids = nan(0,1);
d = nnz(goodidx);
data = nan(0,d);
while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if isempty(s),
    continue;
  end
  ss = regexp(s,'\t','split');
  poscurr = str2double(ss(2:4));
  i = find(pos(:,1) == poscurr(1) & pos(:,2) == poscurr(2) & pos(:,3) == poscurr(3),1);
  if ~isempty(i),
    fprintf('Duplicate x,y,z position found, old blockid: %d, new blockid: %d\n',blockids(i),str2double(ss{1}));
    continue;
  end
  blockids(end+1) = str2double(ss{1});
  pos(end+1,:) = poscurr;
  neuropilnames{end+1} = ss{5};
  neuropilids(end+1) = str2double(ss{6});
  data(end+1,:) = str2double(ss(goodidx));
  
end

fprintf('Read %d rows, %d columns of data\n',size(data));

fclose(fid);