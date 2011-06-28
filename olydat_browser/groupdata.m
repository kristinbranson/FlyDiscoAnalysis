function [gd grpLbl] = groupdata(x,g)
% [gd grpLbl] = groupdata(x,g)
% Group data.
%
% x: data vec. double col vec.
% g: grouping vec. positive int col vec, or col cellstr of same size as x.
% gd: grouped data. see below.
% grpLbl: group labels. int vec or cellstr, corresponding to g.
%
% If g is an int vec: gd is a cell with max(g) number of elements. The ith
% cell of gd is filled with those elements of x where g==i. grpLbl is 1:max(g).
% If g is a cellstr: gd is a cell with numel(unique(g)) number of elements.
% grpLbl is unique(g). The ith cell of gd is filled with those elements of
% x where g==grpLbl{i}.

assert(numel(x)==numel(g));
assert(isnumeric(g) || iscellstr(g));

if isnumeric(g)
  [grpLbl,~,idx] = unique(g);
  Ngrp = numel(grpLbl);
  gd = cell(Ngrp,1);
  for c = 1:Ngrp,
    gd{c} = x(idx==c);
  end
else
    ung = unique(g);
    Ngrp = numel(ung);
    gd = cell(Ngrp,1);
    for c = 1:Ngrp
        gd{c} = x(strcmp(g,ung{c}));
    end
    grpLbl = ung;
end

assert(isequalwithequalnans(sort(cell2mat(gd)),sort(x)));
    
end