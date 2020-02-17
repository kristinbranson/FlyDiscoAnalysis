% function [Mu,Cov,P,Pi,LL]=hmm(X,T,K,cyc,tol);
% 
% Gaussian Observation Hidden Markov Model
%
% X - N x p data matrix
% T - length of each sequence (N must evenly divide by T, default T=N)
% K - number of states (default 2)
% cyc - maximum number of cycles of Baum-Welch (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% Mu - mean vectors
% Cov - output covariance matrix (full, tied across states)
% P - state transition matrix
% Pi - priors
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of Baum-Welch

function [Mu,Cov,P,Pi,LL]=hmm_multiseq(X,K,cyc,tol)

if ~iscell(X),
  X = {X};
end
N = numel(X);
Ts = nan(1,N);
p = size(X{1},2);
for n = 1:N,
  Ts(n) = size(X{n},1);
end

if nargin<4   tol=0.0001; end;
if nargin<3   cyc=100; end;
if nargin<2   K=2; end;

Xmat = cat(1,X{:});
Cov=diag(diag(cov(Xmat)));

Mu=randn(K,p)*sqrtm(Cov)+ones(K,1)*mean(Xmat);

Pi=rand(1,K);
Pi=Pi/sum(Pi);

P=rand(K);
P=rdiv(P,rsum(P));

LL=[];
lik=0;

k1=(2*pi)^(-p/2);

hwait = waitbar(0,'Learning HMM params');

for cycle=1:cyc
  
  %%%% FORWARD-BACKWARD 
  
  Gamma={};
  Gammasum=zeros(1,K);
  %Scale=zeros(T,1);
  %Xi=zeros(T-1,K*K);
  Scale = {};
  Xi = {};
  swait = sprintf('E-Step, iteration %d',cycle);
  if ishandle(hwait),
    waitbar(0,hwait,swait);
  else
    hwait = waitbar(0,swait);
  end
  
  for n=1:N
    
    if ishandle(hwait),
      waitbar(n/(N+1),hwait);
    end
    
    alpha=zeros(Ts(n),K);
    beta=zeros(Ts(n),K);
    gamma=zeros(Ts(n),K);
    B=zeros(Ts(n),K);

    
    iCov=inv(Cov);
    k2=k1/sqrt(det(Cov));
    for i=1:Ts(n)
      for l=1:K
        d=Mu(l,:)-X{n}(i,:);
        B(i,l)=k2*exp(-0.5*d*iCov*d');
      end;
    end;
    
    scale=zeros(Ts(n),1);
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for i=2:Ts(n)
      alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
      scale(i)=sum(alpha(i,:));
      alpha(i,:)=alpha(i,:)/scale(i);
    end;
    
    beta(Ts(n),:)=ones(1,K)/scale(Ts(n));
    for i=Ts(n)-1:-1:1
      beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i); 
    end;
    
    gamma=(alpha.*beta); 
    gamma=rdiv(gamma,rsum(gamma));
    gammasum=sum(gamma);
    
    xi=zeros(Ts(n)-1,K*K);
    for i=1:Ts(n)-1
      t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
      xi(i,:)=t(:)'/sum(t(:));
    end;
    
    %Scale=Scale+log(scale);
    Scale{end+1} = log(scale);
    Gamma{end+1}=gamma;
    Gammasum=Gammasum+gammasum;
    %Xi=Xi+xi;
    Xi{end+1} = xi;
    
  end;
  
  %%%% M STEP 
  
  if ishandle(hwait),
    waitbar(1,hwait,sprintf('\nM-step, iteration %d\n',cycle));
  end

  % outputs
  Mu=zeros(K,p);
  for n = 1:N,
    Mu=Mu + Gamma{n}'*X{n};
  end
  Mu=rdiv(Mu,Gammasum');
  
  % transition matrix 
  sxi = zeros(K,K);
  for n = 1:N,
    sxi=sxi + reshape(rsum(Xi{n}')',[K,K]);
  end
  P=rdiv(sxi,rsum(sxi));
  
  % priors
  Pi=zeros(1,K);
  for i=1:N
    Pi=Pi+Gamma{i}(1,:);
  end
  Pi=Pi/N;
  
  % covariance
  Cov=zeros(p,p);
  for n = 1:N,
    for l=1:K
      d=(X{n}-ones(Ts(n),1)*Mu(l,:));
      Cov=Cov+rprod(d,Gamma{n}(:,l))'*d;
    end
  end;
  Cov=Cov/(sum(Gammasum));
  
  oldlik=lik;
  lik = 0;
  for n = 1:N,
    lik=lik+sum(Scale{n});
  end
  LL=[LL lik];
  fprintf('cycle %i log likelihood = %f ',cycle,lik);  
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    fprintf('violation');
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~isfinite(lik)) 
    fprintf('\n');
    break;
  end;
  fprintf('\n');
end

if ishandle(hwait),
  delete(hwait);
end
