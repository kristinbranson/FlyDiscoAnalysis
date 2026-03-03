function com=border_com(border)

% this calculates the pseudo-'center of mass' of the border.
% the border should be a 2x(n+1) matrix of vertices, with n the number
% of vertices.  the frist row should hold x coords, the second, y coords.
% the (n+1)th col should be identical to the first
% the pseudo-'center of mass' returned is the vector average of the
% vertex coords (ignoring the (n+1)th, of course, since it's a repeat)

n_vertices=size(border,2)-1;
vertices=border(:,1:n_vertices);
com=mean(vertices,2);


