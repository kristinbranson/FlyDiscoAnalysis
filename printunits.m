function s = printunits(units)

if isempty(units.num),
  nums = '1';
else
  nums = units.num{1};
  if numel(units.num) > 1,
    nums = [nums,sprintf('*%s',units.num{2:end})];
  end
end

if isempty(units.den),
  s = nums;
else
  dens = sprintf('/%s',units.den{:});
  s = sprintf('%s%s',nums,dens);
end