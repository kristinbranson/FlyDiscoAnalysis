function [expdirs,movie2flies,fly2movie,nfliespermovie,fly2idx,n2idx] = ...
  SubsetDataStructs(obj,flies,ns)

nflies = length(flies);
nexpdirs = length(ns);
% flyidx i corresponds to fly flies(i)
fly2idx = sparse(ones(1,nflies),flies,1:nflies,1,obj.nflies);
% expidx i corresponds to experiment ns(i)
n2idx = sparse(ones(1,nexpdirs),ns,1:nexpdirs,1,obj.nexpdirs);

expdirs = obj.expdir_bases(ns);
movie2flies = cell(1,nexpdirs);
nfliespermovie = zeros(1,nexpdirs);
for i = 1:length(ns),
  n = ns(i);
  movie2flies{i} = full(fly2idx(obj.movie2flies{n}));
  movie2flies{i}(movie2flies{i} == 0) = [];
  nfliespermovie(i) = length(movie2flies{i});
end
fly2movie = full(n2idx(obj.fly2movie(flies)));
