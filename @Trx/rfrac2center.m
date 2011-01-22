function [x1,y1,x2,y2] = rfrac2center(obj,n,fly)

rfrac = [obj.GetPerFrameData('corfrac_maj',n,fly);obj.GetPerFrameData('corfrac_min',n,fly)];

x_mm = obj.GetPerFrameData('x_mm',n,fly);
y_mm = obj.GetPerFrameData('y_mm',n,fly);
a_mm = obj.GetPerFrameData('a_mm',n,fly);
b_mm = obj.GetPerFrameData('b_mm',n,fly);
theta = obj.GetPerFrameData('theta_mm',n,fly);
x1 = x_mm(1:end-1) + rfrac(1,:).*a_mm(1:end-1).*2.*cos(theta(1:end-1)) - rfrac(2,:).*b_mm(1:end-1).*2.*sin(theta(1:end-1));
y1 = y_mm(1:end-1) + rfrac(1,:).*a_mm(1:end-1).*2.*sin(theta(1:end-1)) + rfrac(2,:).*b_mm(1:end-1).*2.*cos(theta(1:end-1));
x2 = x_mm(2:end) + rfrac(1,:).*a_mm(2:end).*2.*cos(theta(2:end)) - rfrac(2,:).*b_mm(2:end).*2.*sin(theta(2:end));
y2 = y_mm(2:end) + rfrac(1,:).*a_mm(2:end).*2.*sin(theta(2:end)) + rfrac(2,:).*b_mm(2:end).*2.*cos(theta(2:end));


