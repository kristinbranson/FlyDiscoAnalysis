function clr=oppo3_of_x(x)

%  x is a col vector

% normalize x
x=mod(x,1);

% the four cases to be handled
ry= x<0.25 ;
yg= x>=0.25 & x<0.5 ;
gb= x>=0.5 & x<0.75 ;
br= x>=0.75 ;

% the landmark colors
red = [1 0 0] ;
green = [0 0.8 0] ;
blue = [0 0.125 1] ; 
yellow = [1 0.9 0] ;

% the in-between colors
orange=[255 165 0]/255 ;  % red-yellow
purple=[191 0 255]/255 ; % red-blue, "electric purple"
bluegreen = [0, 165, 156]  ;  % Munsell Blue green
chartreuse = [223, 255, 0] ;  % yellow-green

% interpolate between the 'landmark' colors
clr=zeros(length(x),3);
if any(ry~=0)
  y=4*x(ry); r=1-y;
  clr_ry=y*yellow+r*red;  % linear
  clr(ry,:)=clr_ry;
end
if any(yg~=0)
  g=4*x(yg)-1; y=1-g;
  clr_yg=g*green+y*yellow;  % linear
  clr(yg,:)=clr_yg;
end
if any(gb~=0)
  b=4*x(gb)-2; g=1-b;
  % do this one nonlinear so the middle doesn't look dark
  clr_gb=sin(pi/2*b)*blue+sin(pi/2*g)*green;
  clr(gb,:)=clr_gb;
end
if any(br~=0)
  r=4*x(br)-3; b=1-r;
  % do this one nonlinear so the middle doesn't look dark
  clr_br=sin(pi/2*r)*red+sin(pi/2*b)*blue;
  clr(br,:)=clr_br;
end
