function point=nearest_corner(current_point)

% get just the x and y
point=[current_point(1,1) ; current_point(1,2)];

% change to 'corner' coords
point=round(point-0.5)+0.5;
