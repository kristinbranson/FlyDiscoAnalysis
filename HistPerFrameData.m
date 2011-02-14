function [frac,fracless,fracmore,Z,fracframesanalyzed] = HistPerFrameData(data,doanalyze,edges)
    
n = numel(data);
doanalyze = doanalyze & ~isnan(data);
data = data(doanalyze);

Z = nnz(doanalyze);
fracframesanalyzed = Z / n;
    
% histogram with linearly spaced bins
counts = histc(data,edges);
frac = counts / Z;
fracless = frac(1);
fracmore = frac(end-1)+frac(end);
frac = frac(2:end-2);
