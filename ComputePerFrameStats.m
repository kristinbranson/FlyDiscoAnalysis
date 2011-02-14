function [Z,mu,sig,prctiles] = ComputePerFrameStats(data,doanalyze,varargin)

[prctiles_compute] = myparse(varargin,'prctiles_compute',[0,1,5,50,95,99,100]);

doanalyze = doanalyze & ~isnan(data);
data = data(doanalyze);

% number of data points
Z = numel(data);

% mean
mu = mean(data);

% std
sig = std(data,1);

% prctiles
prctiles = prctile(data,prctiles_compute);
