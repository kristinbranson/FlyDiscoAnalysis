function str=alphabetic_label_from_int(i)

% converts ROI ids to alphabetic labels
%  1 ->  a
%  2 ->  b
%    .
%    .
%    .
% 26 ->  z
% 27 -> ba
% 28 -> bb
%    .
%    .
%    .

remnant=i-1;
place=1;
str='';
while true
  digit=mod(remnant,26);
  remnant=floor(remnant/26);
  str(1,place)=digit;
  place=place+1;
  if remnant==0 break; end
end
str=flipdim(str,2);
str=char(str+97);  % want 0->a, 1->b, etc

end
