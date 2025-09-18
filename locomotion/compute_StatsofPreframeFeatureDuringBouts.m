function [datastruct] = compute_StatsofPreframeFeatureDuringBouts(fly,fn,trx,start_indices,end_indices,derivative_flag);
%input fn = name of perframe feature, trx obj, start and end indices of
%bouts in moveie frame of reference

% output = mean value of perframe features during each bout, std of
% perframe features during each bout

dataflycurr = trx.GetPerFrameData(fn,fly);
meanpff = nan(1,numel(start_indices));
stdpff = nan(1,numel(start_indices));
minpff = nan(1,numel(start_indices));
maxpff = nan(1,numel(start_indices));
sumpff = nan(1,numel(start_indices));
n = nan(1,numel(start_indices));

% account for indexing into a derivative perframe feature like velmag
if strcmp('first',derivative_flag)
    end_indices = end_indices-1;
elseif strcmp('second',derivative_flag)
    end_indices = end_indices-2;
end


for i = 1:numel(start_indices)
    if ~isempty(start_indices)
        boutdata = dataflycurr(start_indices(i):end_indices(i));
        meanpff(i) = mean(boutdata);
        stdpff(i) = std(boutdata);
        minpff(i) = min(boutdata);
        maxpff(i) = max(boutdata);
        sumpff(i) = sum(boutdata);
        n(i) = numel(boutdata);


    end

end

datastruct.mean = meanpff;
datastruct.std = stdpff;
datastruct.min = minpff;
datastruct.max = maxpff;
datastruct.sum = sumpff;
datastruct.n = n;
datastruct.pffname = fn;



