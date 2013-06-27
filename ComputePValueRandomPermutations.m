function [fracsmallersamples,fracbiggersamples] = ...
  ComputePValueRandomPermutations(muline,nexpsperset_line,...
  sigmahat_set,sigmahat_exp,nsamples,...
  xexp_control,nexpsperset_control,set2exps_control)

epsilon_compare = .00001;
d = size(xexp_control,1);

nexps_line = sum(nexpsperset_line);
nsets_control = numel(nexpsperset_control);
maxnexps = max(nexpsperset_line);
counts = hist(nexpsperset_line,1:maxnexps);
    
nsmallersamples = zeros(d,1);
nbiggersamples = zeros(d,1);
for samplei = 1:nsamples,

  expidxsample = nan(1,nexps_line);
  n = nan(1,nexps_line);
  isleft = true(1,nsets_control);

  off = 0;
  for minnexps = maxnexps:-1:1,
        
    % choose some sets
    isallowed = nexpsperset_control >= minnexps;
    idxsample = randsample(find(isleft&isallowed),counts(minnexps),false);
    isleft(idxsample) = false;
            
    % choose some experiments
    n(off+1:end+counts(minnexps)*minnexps) = minnexps;
    for seti = 1:counts(minnexps),
      expidxsample(off+1:off+minnexps) = set2exps_control{idxsample(seti)}(randsample(nexpsperset_control(idxsample(seti)),minnexps));
      off = off + minnexps;
    end

  end
        
  badidx = isnan(xexp_control(:,expidxsample));
  expweights = 1./bsxfun(@plus,bsxfun(@times,n,sigmahat_set(:).^2),sigmahat_exp(:).^2);
  expweights = bsxfun(@rdivide,expweights,sum(expweights,2));
  musample = nan(1,d);
  for di = 1:d,
    musample(di) = sum(xexp_control(di,expidxsample(~badidx(di,:))).*expweights(di,~badidx(di,:))) / sum(expweights(di,~badidx(di,:)));
  end
  idx = musample < muline - epsilon_compare;
  nsmallersamples(idx) = nsmallersamples(idx) + 1;
  idx = musample > muline + epsilon_compare;
  nbiggersamples(idx) = nbiggersamples(idx) + 1;
            
end  

fracsmallersamples = nsmallersamples / nsamples;
fracbiggersamples = nbiggersamples / nsamples;