% [ids,mincost] = AssignIdentities_GivenDistributions(x,y,a,b,area,mua,sigmaa,mub,sigmab,muarea,sigmaarea,sigmamotion,vel_dampen)
% assign identities based on area and motion
function [ids,mincost] = AssignIdentities_GivenDistributions(x,y,data,mudata,sigmadata,sigmamotion,vel_dampen,appearanceweight)

isappearanceweight = exist('appearanceweight','var') && ~isempty(appearanceweight);

[nids,T] = size(x);
state2ids = perms(1:nids);
nstates = size(state2ids,1);
ndata = size(data,3);

mudata = reshape(mudata,[nids,1,ndata]);
sigmadata = reshape(sigmadata,[nids,1,ndata]);
constdata = log(sqrt(2*pi)*sigmadata);


% initialize
cost = nan(nstates,nstates,T);
tmpcost = nan(nstates,1);
stateprev = nan(nstates,nstates,T);

% t == 2: nll_appearance for t=1 and t=2 and motion assuming zero velocity
t = 2;
for sc = 1:nstates,
  idsc = state2ids(sc,:);
  % trx(idsc(i),t) corresponds to distribution(i)
  % compute appearance likelihood for current frame
  nll_appearance_c = sum(sum( ((data(idsc,t,:)-mudata)./sigmadata).^2/2 + constdata, 1 ), 3);
      
  if isappearanceweight,
    nll_appearance_c = nll_appearance_c * appearanceweight(t);
  end
  
  for sp = 1:nstates,
    idsp = state2ids(sp,:);
    % trx(idsc(i),t) corresponds to distribution(i)
    % compute appearance likelihood for current frame
    nll_appearance_p = sum(sum( ((data(idsc,t-1,:)-mudata)./sigmadata).^2/2 + constdata, 1 ), 3);
    if isappearanceweight,
      nll_appearance_p = nll_appearance_p * appearanceweight(t-1);
    end

    
    xpred = x(idsp,t-1);
    ypred = y(idsp,t-1);
    nll_motion = sum(( (x(idsc,t) - xpred)./sigmamotion ).^2/2 ...
      + ( (y(idsc,t) - ypred)./sigmamotion ).^2/2);

    cost(sp,sc,t) = nll_motion + nll_appearance_p + nll_appearance_c;
  end
  
end

for t = 3:T,
  
  for sc = 1:nstates,
    
    idsc = state2ids(sc,:);
    % trx(idsc(i),t) corresponds to distribution(i)
    
    % compute appearance likelihood for current frame
    nll_appearance_c = sum(sum( ((data(idsc,t,:)-mudata)./sigmadata).^2/2 + constdata, 1 ), 3);
    
    if isappearanceweight,
      nll_appearance_c = nll_appearance_c * appearanceweight(t);
    end
    
    for sp = 1:nstates,
      
      idsp = state2ids(sp,:);
      % trx(idsp(i),t-1)
      
      for spp = 1:nstates,
        idspp = state2ids(spp,:);
        % trx(idspp(i),t-2)
        
        xpred = (2-vel_dampen)*x(idsp,t-1) - (1-vel_dampen)*x(idspp,t-2);
        ypred = (2-vel_dampen)*y(idsp,t-1) - (1-vel_dampen)*y(idspp,t-2);
        nll_motion = sum(( (x(idsc,t) - xpred)./sigmamotion ).^2/2 ...
          + ( (y(idsc,t) - ypred)./sigmamotion ).^2/2);
        
        tmpcost(spp) = nll_motion + cost(spp,sp,t-1);
        
      end
      
      [mintmpcost,best_spp] = min(tmpcost);
      cost(sp,sc,t) = mintmpcost + nll_appearance_c;
      stateprev(sp,sc,t) = best_spp;
      
    end
  end
  
end

cost = reshape(cost,[nstates*nstates,T]);
%stateprev = reshape(stateprev,[nstates*nstates,T]);

% choose the best last two states
[mincost,s] = min(cost(:,T),[],1);
[sp,sc] = ind2sub([nstates,nstates],s);
ids = nan(nids,T);
ids(:,T) = state2ids(sc,:);
ids(:,T-1) = state2ids(sp,:);
spp = stateprev(sp,sc,T);
ids(:,T-2) = state2ids(spp,:);

for t = T-1:-1:3,
  sc = sp;
  sp = spp;
  spp = stateprev(sp,sc,t);
  ids(:,t-2) = state2ids(spp,:);
end
