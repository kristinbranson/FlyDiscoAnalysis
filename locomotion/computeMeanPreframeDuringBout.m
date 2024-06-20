function [meanperframe,stdperframe] = computeMeanPreframeDuringBout(fly,fn,trx,start_indices,end_indices);
%input fn = name of perframe feature, trx obj, start and end indices of
%bouts in moveie frame of reference

% output = mean value of perframe features during each bout, std of
% perframe features during each bout

dataflycurr = trx.GetPerFrameData(fn,fly);
offcurrfly = trx(fly).off;
for i = 1:numel(start_indices)
    meanperframe(i) = mean(dataflycurr(start_indices(i)+offcurrfly:end_indices(i)+offcurrfly));
    stdperframe(i) = std(dataflycurr(start_indices(i)+offcurrfly:end_indices(i)+offcurrfly));


end
