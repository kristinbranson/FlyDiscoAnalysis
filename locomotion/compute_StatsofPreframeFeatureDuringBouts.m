function [datastruct] = compute_StatsofPreframeFeatureDuringBouts(fly,fn,trx,start_indices,end_indices,dataflycurr)
%input fn = name of perframe feature, trx obj, start and end indices of
%bouts in moveie frame of reference
% optional: dataflycurr = pre-fetched perframe data to avoid repeated trx.GetPerFrameData calls

% output = mean value of perframe features during each bout, std of
% perframe features during each bout

if nargin < 6 || isempty(dataflycurr)
    dataflycurr = trx.GetPerFrameData(fn,fly);
end
% deal with empty case - return NaN OK?
if isempty(start_indices)
    start_indices = 1;    
    meanpff = nan(1,numel(start_indices));
    stdpff = nan(1,numel(start_indices));
    minpff = nan(1,numel(start_indices));
    maxpff = nan(1,numel(start_indices));
    sumpff = nan(1,numel(start_indices));
    n = nan(1,numel(start_indices));
    start_indices = [];
    allboutdata = [];
else
    meanpff = nan(1,numel(start_indices));
    stdpff = nan(1,numel(start_indices));
    minpff = nan(1,numel(start_indices));
    maxpff = nan(1,numel(start_indices));
    sumpff = nan(1,numel(start_indices));
    n = nan(1,numel(start_indices));

    % end_indices are exclusive (first frame after the bout). Perframe
    % features are aligned to frame i (raw: value at i; first derivative:
    % forward diff i->i+1), so the last in-bout index is end-1. Clamp to the
    % available data for the terminal bout, where the final frame/transition
    % has no entry (e.g. a length n-1 derivative on the last frame).
    end_indices = min(end_indices - 1, numel(dataflycurr));

    allboutdata = {};%
    for i = 1:numel(start_indices)
        boutdata = dataflycurr(start_indices(i):end_indices(i));
        allboutdata{i} = boutdata;%
        meanpff(i) = mean(boutdata);
        stdpff(i) = std(boutdata);
        minpff(i) = min(boutdata);
        maxpff(i) = max(boutdata);
        sumpff(i) = sum(boutdata);
        n(i) = numel(boutdata);
    end
    allboutdata = horzcat(allboutdata{:}); % AR 2026019 fixing boutdata only saving data for last bout

end
datastruct.data = allboutdata; 
datastruct.mean = meanpff;
datastruct.std = stdpff;
datastruct.min = minpff;
datastruct.max = maxpff;
datastruct.sum = sumpff;
datastruct.n = n;
datastruct.pffname = fn;


