function result = outsideArenaMaskFromRegistrationData(registration_data, nc, nr)

xc = registration_data.circleCenterX ;
yc = registration_data.circleCenterY ;
r = registration_data.circleRadius ;
chamber_count = numel(xc) ;
xc3 = reshape(xc,[1 1 chamber_count]) ;
yc3 = reshape(yc,[1 1 chamber_count]) ;

[xgrid, ygrid] = meshgrid(1:nc, 1:nr);
isOutsideChamber = ((xgrid-xc3).^2 + (ygrid-yc3).^2) >= r.^2 ;
result = all(isOutsideChamber,3) ;
