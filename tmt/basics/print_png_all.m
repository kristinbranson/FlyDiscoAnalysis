function f(hs, res)

if nargin<1 || isempty(hs)
  hs=get(0,'Children')';
end
if nargin<2 
    res = [] ;
end
% print em
for i=hs
  print_png(i, [], res);
  if i~=hs(end) fprintf(1,'\n'); end
end
