function result = smooth_colormap_in_six_spans(n_colors, cmap_fn)

% generate a non-smooth colormap on a very fine grid
% calculate the path length, and from that get the phase for each
% point on the fine grid
n_samples=6000;  % want divisible by 6
x=linspace(0,1,n_samples+1)';
clr=cmap_fn(x);
clr_lab=srgb2lab(clr);
ds=dist_lab(clr_lab(1:end-1,:),clr_lab(2:end,:));
s=[0 ; cumsum(ds)];
%phase=s/s(end);  % normalized path length == phase in cycles

% normalize in a special way, to leave certain points where we want
% them
s1=interp1(x,s,1/6,'linear');
s2=interp1(x,s,2/6,'linear');
s3=interp1(x,s,3/6,'linear');
s4=interp1(x,s,4/6,'linear');
s5=interp1(x,s,5/6,'linear');
s6=interp1(x,s,6/6,'linear');
span1= x<1/6;
span2= x>=1/6 & x<2/6;
span3= x>=2/6 & x<3/6;
span4= x>=3/6 & x<4/6;
span5= x>=4/6 & x<5/6;
span6= x>=5/6;
phase=zeros(size(s));
phase(span1)=1/6/s1*s(span1);
phase(span2)=1/6/(s2-s1)*(s(span2)-s1)+1/6;
phase(span3)=1/6/(s3-s2)*(s(span3)-s2)+2/6;
phase(span4)=1/6/(s4-s3)*(s(span4)-s3)+3/6;
phase(span5)=1/6/(s5-s4)*(s(span5)-s4)+4/6;
phase(span6)=1/6/(s6-s5)*(s(span6)-s5)+5/6;

% make a colormap with inter-color spacings equal to circum/n_colors
phase_samples=linspace(0,1,n_colors+1)';
phase_samples=phase_samples(1:end-1);
x_samples=interp1(phase,x,phase_samples,'linear');
result=cmap_fn(x_samples);

