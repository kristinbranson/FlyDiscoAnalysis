function [sety,expy] = CollectSetDataGivenGenotype(fn,genotype,setstats,expstats,metadata)

setidx = find(strcmp({setstats.metadata.genotype},genotype));
nsetscurr = numel(setidx);

expy = cell(1,nsetscurr);
sety = nan(1,nsetscurr);
for j = 1:nsetscurr,
  seti = setidx(j);
  setname = setstats.metadata(seti).set;
  expidx = strcmp({metadata.set},setname);
  
  expy{j} = sort(expstats.(fn)(expidx));
  expy{j}(isnan(expy{j})) = [];
  sety(j) = setstats.means.(fn)(seti);
end

badidx = isnan(sety);
sety(badidx) = [];
expy(badidx) = [];
