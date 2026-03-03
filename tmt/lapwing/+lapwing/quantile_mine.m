function x_p=quantile_mine(x,p)

% Calculate quantiles of the data in x.  (E.g. the 0.5 quantile is the
% median, the 0.75 quantile is the 75th percentile.)
% x a col, p a col.  x_p same size as p

% We use a very stupid algorithm, because the alternatives are worse.

x_sorted=sort(x);
n=length(x);
P=linspace(0,1,n)';
x_p=interp1q(P,x_sorted,p);

end
