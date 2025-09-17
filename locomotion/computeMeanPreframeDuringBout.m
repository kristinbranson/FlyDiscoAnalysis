function [meanperframe,stdperframe] = computeMeanPreframeDuringBout(fly,fn,trx,start_indices,end_indices,derivative_flag);
%input fn = name of perframe feature, trx obj, start and end indices of
%bouts in moveie frame of reference

% output = mean value of perframe features during each bout, std of
% perframe features during each bout

dataflycurr = trx.GetPerFrameData(fn,fly);
meanperframe = nan(1,numel(start_indices));
stdperframe = nan(1,numel(start_indices));

% account for indexing into a derivative perframe feature like velmag
if strcmp('first',derivative_flag)
    end_indices = end_indices-1;
elseif strcmp('second',derivative_flag)
    end_indices = end_indices-2;
end


for i = 1:numel(start_indices)
    if ~isempty(start_indices)
        meanperframe(i) = mean(dataflycurr(start_indices(i):end_indices(i)));
        stdperframe(i) = std(dataflycurr(start_indices(i):end_indices(i)));
    end

end
