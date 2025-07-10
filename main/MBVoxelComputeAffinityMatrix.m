function [W,w] = MBVoxelComputeAffinityMatrix(data,D,ndimeff,cproplaplace,sigma)

ndata = size(data,1);

% find edges
idx = ~isinf(D);
nedges = nnz(idx);
[r,c] = find(~isinf(D));
dloc = round(double(D(idx)));
maxdloc = max(dloc);

% what should be the b for each edge?
bfactor = ones(size(dloc));
bfactorcurr = 1;
for i = 1:maxdloc,
  bfactor(dloc == i) = bfactorcurr;
  bfactorcurr = bfactorcurr*(1+cproplaplace/i);
end
dint = sum(abs(data(r,:)-data(c,:)),2);

if numel(sigma) == 1,
  w = exp(-dint./bfactor/sigma)./bfactor.^ndimeff;
else
  sigmaedge = sqrt(sigma(r).*sigma(c));
  w = exp(-dint./bfactor./sigmaedge)./bfactor.^ndimeff;
end
  
W = sparse(r,c,w,ndata,ndata,nedges);