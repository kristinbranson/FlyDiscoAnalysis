sig = 5;
x1 = randn(10,10000)*sig;
x2 = randn(10,10000)*sig + x1;
x3 = randn(10,10000)*sig + x2;

cov(x3')

%%

mus = 1:100;
ns = 1:15;

b = nan(numel(mus),numel(ns));

for mui = 1:numel(mus),
  
  mu = mus(mui);
    
  x11 = exprnd(mu,1,100000);
  x12 = exprnd(mu,1,100000);
  x1 = x11-x12;
  
  b(mui,1) = sqrt(var(x1,1,2)/2);

  for ni = 1:numel(ns),
  
    n = ns(ni);
    x2 = x1;
  
    for i = 1:n,
      x21 = exprnd(mu,1,100000);
      x22 = exprnd(mu,1,100000);
      x2 = x21-x22 + x2;
    end
    
    b(mui,ni+1) = sqrt(var(x2,1,2)/2);
    
  end
end

bratio = b(:,2:end)./b(:,1:end-1);

%%
clf;
%plot([0,ns],b,'.-');
plot(ns,bratio,'.-');
legend(num2str(mus'));

hold on;
% deg = 1;
% p = polyfit(repmat(1./ns,[numel(mus),1]),bratio,deg);
% y = p(end);
% for i = 1:deg,
%   y = y + ns.^(-(deg-i+1))*p(i);
% end

c = sum(sum(bratio-1))./(numel(mus)*sum(1./ns));

y = c./ns + 1;
h = plot(ns,y,'k--');
