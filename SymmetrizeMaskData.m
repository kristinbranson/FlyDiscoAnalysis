function maskdata = SymmetrizeMaskData(maskdata)

if isfield(maskdata,'leg_symmetric') && isfield(maskdata,'mask_symmetric'),
  return;
end

% symmetrize
fprintf('Creating symmetric version of the mask...\n');
mleft = regexp(maskdata.leg,'^(.*)_L$','once','tokens');
isleft = ~cellfun(@isempty,mleft);
idxleft = find(isleft);
mleft = cellfun(@(x) x{1},mleft(isleft),'UniformOutput',false);
mright = regexp(maskdata.leg,'^(.*)_R$','once','tokens');
isright = ~cellfun(@isempty,mright);
idxright = find(isright);
mright = cellfun(@(x) x{1},mright(isright),'UniformOutput',false);
isneither = ~isleft & ~isright;
idxneither = find(isneither);
maskdata.mask_symmetric = zeros(size(maskdata.mask),'uint8');
ncompartments = numel(idxleft)+nnz(isneither);
maskdata.leg_symmetric = cell(1,ncompartments);
[~,idxleft2right] = ismember(mleft,mright);
compi = 1;
for i = 1:numel(idxneither),
  maskdata.mask_symmetric(maskdata.mask == idxneither(i)) = compi;
  maskdata.leg_symmetric{compi} = maskdata.leg{idxneither(i)};
  compi = compi + 1;
end
for i = 1:numel(idxleft),
  ileft = idxleft(i);
  iright = idxright(idxleft2right(i));
  maskdata.mask_symmetric(maskdata.mask == ileft | maskdata.mask == iright) = compi;
  maskdata.leg_symmetric{compi} = mleft{i};
  compi = compi + 1;
end
maskdata.mask_symmetric = permute(maskdata.mask_symmetric,[2,1,3]);

fprintf('Computing spatial extent of each compartment...\n');
% get the spatial extent of each compartment
maskdata.xlims = nan(ncompartments,2);
maskdata.ylims = nan(ncompartments,2);
maskdata.zlims = nan(ncompartments,2);
maskdata.npx = nan(1,ncompartments);
for i = 1:ncompartments,
  is3 = any(any(maskdata.mask_symmetric==i,1),2);
  is2 = any(any(maskdata.mask_symmetric==i,1),3);
  is1 = any(any(maskdata.mask_symmetric==i,2),3);
  maskdata.xlims(i,:) = [find(is2,1),find(is2,1,'last')];
  maskdata.ylims(i,:) = [find(is1,1),find(is1,1,'last')];
  maskdata.zlims(i,:) = [find(is3,1),find(is3,1,'last')];
  maskdata.npx(i) = nnz(maskdata.mask_symmetric==i);
end


% get the x-y outlines of each mask
fprintf('Computing boundaries of each compartment...\n');
maskdata.boundaries = struct;
for i = 1:ncompartments,
  tmp = any(maskdata.mask_symmetric==i,3);
  b = bwboundaries(tmp);
  maskdata.boundaries(i).y = [];
  maskdata.boundaries(i).x = [];
  for j = 1:numel(b),
    maskdata.boundaries(i).y = [maskdata.boundaries(i).y;b{j}(:,1);nan];
    maskdata.boundaries(i).x = [maskdata.boundaries(i).x;b{j}(:,2);nan];
  end
  maskdata.boundaries(i).y(end) = [];
  maskdata.boundaries(i).x(end) = [];
end