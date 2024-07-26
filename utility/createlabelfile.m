% create label file for behavior 
% behaviorlabel = 'jump';
% array = isjumping;
% 
% offset = trx.firstframes(i);
% 
% fly = i; 
function [flies, t0s, t1s, names, off] = createlabelfile(behaviorlabel,array,offset, fly)


transitions = diff(array);

starts = find(transitions == 1);
starts = starts+1+offset;
if array(1) == 1
    t0s(1) = 1+offset;
    t0s = [t0s, starts];
else
    t0s = starts;
end

ends = find(transitions == -1);
ends = ends;
t1s = ends+offset;
if array(end) == 1
    t1s = [t1s, length(array)+offset];
end
% t0s
% t1s

flies = fly; 
names = {behaviorlabel};
names = repmat(names,1,length(t0s));
off = 1-offset;
%timestamp = now;

end




