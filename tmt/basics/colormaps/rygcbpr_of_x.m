function clr = rygcbpr_of_x(x)

% red-yellow-green-cyan-blue-purple-red

%  x is a col vector

% normalize x
x=mod(x,1);

% the four cases to be handled
ry= x<1/6 ;
yg= x>=1/6 & x<2/6 ;
gc= x>=2/6 & x<3/6 ;
cb= x>=3/6 & x<4/6 ;
bp= x>=4/6 & x<5/6 ;
pr= x>=5/6 ;

% the landmark colors
red = [1 0 0] ;
yellow = [1 0.9 0] ;
green = [0 0.8 0] ;
cyan = [0 0.85 0.85]  ;
blue = [0 0.125 1] ; 
purple=[0.75 0 1] ; % "electric purple"

% % the in-between colors
% orange=[1 0.65 0] ;  % red-yellow
% chartreuse = [0.87 1 0] ;  % yellow-green

% interpolate between the 'landmark' colors
clr=zeros(length(x),3);
if any(ry~=0)
  y=6*x(ry); r=1-y;
  clr_ry=y*yellow+r*red;  % linear
  clr(ry,:)=clr_ry;
end
if any(yg~=0)
  g=6*x(yg)-1; y=1-g;
  clr_yg=g*green+y*yellow;  % linear
  clr(yg,:)=clr_yg;
end
if any(gc~=0)
  c=6*x(gc)-2; g=1-c;
  clr_gc=g*green+c*cyan;  % linear
  clr(gc,:)=clr_gc;
end
if any(cb~=0)
  b=6*x(cb)-3; c=1-b;
  clr_cb=b*blue+c*cyan;  % linear
  clr(cb,:)=clr_cb; 
end
if any(bp~=0)
  p=6*x(bp)-4; b=1-p;
  clr_bp=b*blue+p*purple;  % linear
  clr(bp,:)=clr_bp; 
end
if any(pr~=0)
  r=6*x(pr)-5; p=1-r;
  clr_pr=r*red+p*purple;  % linear
  clr(pr,:)=clr_pr; 
end

