% states = first_order_viterbi(x,nstates,log_x_trans_fun,log_s_trans_fun)
% 
% Dynamic programming solution for finding the highest probability
% path. 
%
% Assumed Model:
%
% p(x_t | x_{t-1}, s_t) = x_trans_fun(x_t,x_{t-1},s_t)
% P(s_t | s_{t-1}) = s_trans_fun(s_t,s_{t-1})
%
% \prod_{\tau=1}^t p(x_\tau | x_{\tau-1}, s_{\tau}) P(s_\tau | s_{\tau-1})
% = \prod_{\tau=1}^t p(x_\tau, s_\tau | x_{\tau-1}, s_{\tau-1})
% = p(x_{1:t},s_{1:t})
%
% Recursion:
%
% let ll(t,a) = maximum log likelihood over all possible previous
% states s_{1:t-1} of x_{1:t}, s_t = a:
% ll(t,a) = max_{s_{1:t-1}} log p(x_{1:t}, s_{1:t-1}, s_t = a)
% = max_{s_{t-1}} max_{s_{1:t-2}} [p(x_{1:t-1},s_{1:t-1}) * 
%                                  p(x_t | x_{t-1}, s_t=a) *
%                                  P(s_t=a | s_{t-1})]
% = max_{s_{t-1}} [ p(x_t | x_{t-1}, s_t=a) * P(s_t=a | s_{t-1}) * 
%                   max_{s_{1:t-2}} [p(x_{1:t-1},s_{1:t-1})] ]
% = max_{s_{t-1}}[p(x_t|x_{t-1},s_t=a) * P(s_t=a|s_{t-1}) * ll(t-1,s_{t-1})]

function [states,normscore] = first_order_viterbi(x,nstates,log_x_trans_fun,...
  log_s_trans_fun)

T = size(x,1);

% initialize DP
% score
ll = zeros(T,nstates);
% store the path by storing the previous state
prev = nan(T,nstates);

% ll(1,a) is initialized to 1/nstates
ll(1,:) = 1/nstates;

% compute iteratively
for t = 2:T,

  % get the score for each of the possibel current states
  scorex = zeros(nstates,1);
  score = zeros(nstates,1);
  for a = 1:nstates,
    scorex(a) = log_x_trans_fun(x(t,:),x(t-1,:),a);
  end
  %if all(scorex < eps),
  %    scorex(:) = 1;
  %end
  %scorex = scorex / sum(scorex);
  %badidx = scorex < .000001;
  %scorex(badidx) = .000001;
  %scorex = log(scorex);
  %scorex(badidx) = -inf;
  
  for a = 1:nstates,
    for b = 1:nstates,      
      score(b) = ll(t-1,b) + scorex(a) + log_s_trans_fun(a,b);
    end
    [ll(t,a),prev(t,a)] = max(score);
  end
  
end

% get the best chain
states = nan(T,1);
% choose the best last state
[score,states(T)] = max(ll(T,:));
for t = T-1:-1:1,
  states(t) = prev(t+1,states(t+1));
end
normscore = score / T;