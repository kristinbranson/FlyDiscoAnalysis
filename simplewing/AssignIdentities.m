function [idscurr,mudata,sigmadata,sigmamotion,costcurr,iter] = ...
  AssignIdentities(x,y,data,varargin)

[vel_dampen,sigmamotion,niters,appearanceweight,sortbyidx] = myparse(varargin,'vel_dampen',0,...
  'sigmamotion',[],'niters',100,'appearanceweight',[],'sortbyidx',1);
estimatesigmamotion = isempty(sigmamotion);

[nids,T] = size(x);
ndata = size(data,3);

isappearanceweight = ~isempty(appearanceweight);

% initialize distribution of sizes by setting id1 to be the smaller in
% every frame, id2 to be the larger
datasort = nan(size(data));
[~,order] = sort(data(:,:,sortbyidx),1);
for i = 1:ndata,
  datasort(:,:,i) = data(sub2ind([nids,T,ndata],order,repmat(1:T,[nids,1]),repmat(i,[nids,T])));
end
mudata = nan(nids,ndata);
sigmadata = nan(nids,ndata);

if isappearanceweight
  
  zappearance = sum(appearanceweight);

  for i = 1:ndata,
    mudata(:,i) = sum(bsxfun(@times,datasort(:,:,i),appearanceweight),2)/zappearance;
    sigmadata(:,i) = sqrt(sum(bsxfun(@times,bsxfun(@minus,datasort(:,:,i),mudata(:,i)).^2,appearanceweight),2)/zappearance);
  end
  
else
  mudata = reshape(mean(dataasort,2),[nids,ndata]);
  sigmadata = reshape(std(datasort,1,2),[nids,ndata]);
end

sigmadata = max(sigmadata,.0001);

state2ids = perms(1:nids);
nstates = size(state2ids,1);

if any(isnan(mudata(:))) || any(isnan(sigmadata(:))),
  error('NaNs found in initialization');
end

if estimatesigmamotion,

  % initialize distribution of sigmamotion by taking the min motion per
  % frame
  err = inf(1,T-1);
  for sp = 1:nstates,
    idsp = state2ids(sp,:);
    for spp = 1:nstates,
      idspp = state2ids(spp,:);
      xpred = [x(idsp,1),(2-vel_dampen)*x(idsp,2:end-1) - (1-vel_dampen)*x(idspp,1:end-2)];
      ypred = [y(idsp,1),(2-vel_dampen)*y(idsp,2:end-1) - (1-vel_dampen)*y(idspp,1:end-2)];
      errcurr = sum((x(:,2:end)-xpred).^2 + (y(:,2:end)-ypred).^2,1);
      err = min(err,errcurr);
    end
  end
  
  sigmamotion = sqrt( sum(err)/(2*nids*T) );
  
  if isnan(sigmamotion),
    error('NaN found for sigmamotion');
  end
  
end

% 
% state2ids = perms(1:nids);
% nstates = size(state2ids,1);

for iter = 1:niters,
 
  [idscurr,costcurr] = AssignIdentities_GivenDistributions(x,y,data,mudata,sigmadata,sigmamotion,vel_dampen,appearanceweight);

  if iter > 1 && all(idsprev(:) == idscurr(:)),
    break;
  end
  
  if iter > 1 && costcurr > costprev,
    warning('Increase in cost at iteration %d by %f%%',(costcurr-costprev)/costprev);
  end
  
  % update distributions
  xsort = x(sub2ind([nids,T],idscurr,repmat(1:T,[nids,1])));
  ysort = y(sub2ind([nids,T],idscurr,repmat(1:T,[nids,1])));
  for i = 1:ndata,
    datasort(:,:,i) = data(sub2ind([nids,T,ndata],idscurr,repmat(1:T,[nids,1]),repmat(i,[nids,T])));
  end
  
  if isappearanceweight
  
    for i = 1:ndata,
      mudata(:,i) = sum(bsxfun(@times,datasort(:,:,i),appearanceweight),2)/zappearance;
      sigmadata(:,i) = sqrt(sum(bsxfun(@times,bsxfun(@minus,datasort(:,:,i),mudata(:,i)).^2,appearanceweight),2)/zappearance);
    end
    
  else
    mudata = reshape(mean(dataasort,2),[nids,ndata]);
    sigmadata = reshape(std(datasort,1,2),[nids,ndata]);
  end
  
  if any(isnan(mudata(:))),
    error('NaNs found at iteration %d',iter);
  end
  
  if estimatesigmamotion,
    xpred = [xsort(:,1),(2-vel_dampen)*xsort(:,2:end-1) - (1-vel_dampen)*xsort(:,1:end-2)];
    ypred = [ysort(:,1),(2-vel_dampen)*ysort(:,2:end-1) - (1-vel_dampen)*ysort(:,1:end-2)];  
    sigmamotion = sqrt( sum(sum((xsort(:,2:end)-xpred).^2 + (ysort(:,2:end)-ypred).^2,1),2) / (2*nids*T) );
    
    if isnan(sigmamotion),
      error('NaN found for sigmamotion at iteration %d',iter);
    end 
  end
  
  if isnan(costcurr),
    error('NaN found for cost at iteration %d',iter);
  end
    
  idsprev = idscurr;
  costprev = costcurr;
  
end

% make sure that mudata is still in the right order
[~,order] = sort(mudata(:,sortbyidx));
mudata = mudata(order,:);
sigmadata = sigmadata(order,:);
idscurr = idscurr(order,:);

datasort = nan(size(data));
[~,order] = sort(data(:,:,sortbyidx),1);
for i = 1:ndata,
  datasort(:,:,i) = data(sub2ind([nids,T,ndata],order,repmat(1:T,[nids,1]),repmat(i,[nids,T])));
end