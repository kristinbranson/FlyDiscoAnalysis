function playufmf(filename)
if ~exist('filename', 'var') || isempty(filename) ,
  playfmf() ;
else
  playfmf([],filename) ;
end
