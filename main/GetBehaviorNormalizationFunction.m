% fcn = GetBehaviorNormalizationFunction( X, y, ...)
% X is N x D
% y is N x 1
function [fcn,coeffs] = GetBehaviorNormalizationFunction(X, y, coeffs )

if nargin < 3,
  
  [N,D] = size(X);
  y = y(:);
  if numel(y) ~= N,
    error('Sizes of X and y do not match');
  end
  
  % D+1 x 1
  coeffs = regress(y,[ones(N,1),X]);
  
end

fcn = @(yy,XX) yy - max(0,[ones(size(XX,1),1),XX]*coeffs);

%fcn = @(yy,XX) BehaviorNormalizationFunction(yy,XX,coeffs);

% function z = BehaviorNormalizationFunction(y,X,coeffs)
% 
% [N,D] = size(X);
% y = y(:);
% if numel(y) ~= N,
%   error('Sizes of X and y do not match');
% end
% 
% ypred = [ones(N,1),X]*coeffs;
% z = y - ypred;