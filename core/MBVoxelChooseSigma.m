function [sigma,sigmalocal] = MBVoxelChooseSigma(D,data,kneighbors)

ndata = size(data,1);

idx = ~isinf(D);
[r,c] = find(idx);
dint = sum(abs(data(r,:)-data(c,:)),2);
dloc = double(D(idx));

% choose sigma
% compute distance to kneighbors-nearest neighbor
sigmalocal = nan(ndata,1);
for i = 1:ndata,
  dneighbors = sort(dint(r==i & round(dloc)<=1));
  if isempty(dneighbors),
    sigmalocal(i) = min(dint(r==i));
    if isinf(sigmalocal(i)),
      sigmalocal(i) = nan;
    end
    continue;
  end
  if numel(dneighbors) < kneighbors,
    sigmalocal(i) = dneighbors(end);
  else
    sigmalocal(i) = dneighbors(kneighbors);
  end
end
sigma = nanmedian(sigmalocal);
