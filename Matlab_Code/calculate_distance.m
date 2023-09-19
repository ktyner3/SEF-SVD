function [distance] = calculate_distance(x,y)
% Calculate the distance between a dipole point (x) in 3D space with other
% points (y) in 3D space.
% x should be a 1x3 vector, and y should be an nx3 matrix.
% Positions should be in x,y,z coordinate system.
% Function will calculate distance of point x to each point in y.
[m,~] = size(y);
distance = zeros(m,1);
xx = x(1,1); % x in x direction
xy = x(1,2); % x in y direction
xz = x(1,3); % x in z direction

for ii = 1:m
    yx = y(ii,1); % y in x direction
    yy = y(ii,2); % y in y direction
    yz = y(ii,3); % y in z direction

    pt_dist = sqrt( ((xx-yx)^2) + ((xy-yy)^2) + ((xz-yz)^2) );

    distance(ii,1) = pt_dist;
end

end
