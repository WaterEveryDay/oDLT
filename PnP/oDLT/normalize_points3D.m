function [Xtilde, T, T_inv] = normalize_points3D(X)
    % INPUT:
    % X: nx3 matrix, where each row is a 3D point [x_i, y_i, z_i]
    %
    % OUTPUT:
    % Xtilde: nx3 matrix, normalized points
    % T: 4x4 similarity transformation matrix

    % Step 1: Compute the centroid of the points
    centroid = mean(X, 1);  % [c_x, c_y, c_z]

    % Step 2: Translate the points so the centroid is at the origin
    X_translated = X - centroid;  % Subtract centroid from each point

    % Step 3: Compute the average squared distance of the translated
    % points from the origin
    d2_avg = mean(sum(X_translated(:,1:3).^2, 2));

    % Step 4: Scale the points so that the average distance is sqrt(2)
    scale_factor = sqrt(3 / d2_avg);
    Xtilde = scale_factor * X_translated;

    % Step 5: Construct the similarity transformation matrix T
    % T includes translation and scaling
    T = [scale_factor, 0, 0, -scale_factor * centroid(1);
         0, scale_factor, 0, -scale_factor * centroid(2);
         0, 0, scale_factor, -scale_factor * centroid(3);
         0, 0, 0, 1];

    T_inv = [1/scale_factor, 0, 0, centroid(1);
             0, 1/scale_factor, 0, centroid(2);
             0, 0, 1/scale_factor, centroid(3);
             0, 0, 0, 1];

end