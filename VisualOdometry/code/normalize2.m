function [T,normpoints] = normalize2(points)
%% the function normalizes the points by translating and scalling those points in a way that the centroid of all the points be at origin and the mean distance of all the points to the origin is 2^(1/2)

mean_points = mean(points);
tx = mean_points(1,1);
ty = mean_points(1,2);

newpoints = points - repmat(mean_points,size(points, 1),1);
    temp = sqrt((mean(sum((newpoints.^2),2))));
    scale = sqrt(2)/temp;
    T = [scale 0 0 ;0 scale 0; 0 0 1] * [1 0 -tx; 0 1 -ty; 0 0 1];
    
    normpoints = zeros(size(points));
     for i = 1: size(points, 1)
        tempp = points(i, :);
        XY = T * [tempp(1, 1); tempp(1, 2); 1 ];

        normpoints(i, :) = [ XY(1) XY(2)];
     end
     normpoints = normpoints;
end

