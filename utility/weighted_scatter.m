% [mu,S] = weighted_scatter(x,mu,w)
% x is n x d, w is n x 1
function S = weighted_scatter(x,mu,w)

%[n,d] = size(x);

z = sum(w);

diffs = bsxfun(@minus,x,mu);
diffs = bsxfun(@times,diffs,sqrt(w));
S = (diffs'*diffs)/z;
