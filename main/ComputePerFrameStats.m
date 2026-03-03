function [Z,mu,sig,prctiles] = ComputePerFrameStats(data,doanalyze,varargin)

[prctiles_compute] = myparse(varargin,'prctiles_compute',[0,1,5,50,95,99,100]);

doanalyze = doanalyze & ~isnan(data) & ~isinf(data);
data = data(doanalyze);

% number of data points
Z = numel(data);

% mean
mu = mean(data);

% std
sig = std(data,1);

% prctiles
if isempty(prctiles_compute),
  prctiles = nan(0,1);
else
  prctiles = prctile(data,prctiles_compute);
end
