function x = padgrab(trx,fn,flies,f0,f1,padv)

if nargin < 6,
  padv = nan;
end
x = trx(flies).(fn);
for i = 1:numel(flies),
  x{i} = padgrab(x{i},padv,1,size(x{i},1),f0,f1);
end