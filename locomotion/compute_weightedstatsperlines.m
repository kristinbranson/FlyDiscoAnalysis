function [mean_weighted,SE_weighted] = compute_weightedstatsperlines(moviemeans,moviestds,Ns)

% input = 
% meanvalues = (nexps,bins)
% Ns = (nexps,bins)
mean_weighted = nan(1,size(moviemeans,1));
SE_weighted = nan(1,size(moviemeans,1));


mean_weighted = sum(moviemeans.*Ns,1,'omitnan')./sum(Ns,1,'omitnan');


sigma_weighted = sum(Ns.*(moviestds.^2 + (moviemeans-mean_weighted).^2),1,'omitnan')./sum(Ns,1,'omitnan');

SE_weighted = sqrt(sigma_weighted)./sqrt(sum(Ns,1,'omitnan'));
    