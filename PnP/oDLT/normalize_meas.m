function [Utilde, T, T_inv] = normalize_meas(U)
    % INPUT:
    % U: nx3 or nx2 matrix, where each column is a 2D point [x_i, y_i]
    %
    % OUTPUT:
    % Utilde: nx2 matrix, normalized points
    % T: 3x3 similarity transformation matrix

    % Step 1: Compute the centroid of the points
    centroid = mean(U, 1);  % [c_x, c_y]

    % Step 2: Translate the points so the centroid is at the origin
    U_translated = U - centroid;  % Subtract centroid from each point

    % Step 3: Compute the average distance of the translated points from the origin
    d2_avg = mean(sum(U_translated(:,1:2).^2, 2));

    % Step 4: Scale the points so that the average distance is sqrt(2)
    scale_factor = sqrt(2 / d2_avg);
    Utilde = scale_factor * U_translated;
    Utilde(:,3) = ones(size(U, 1), 1);

    % Step 5: Construct the similarity transformation matrix T
    % T includes translation and scaling
    T = [scale_factor, 0, -scale_factor * centroid(1);
         0, scale_factor, -scale_factor * centroid(2);
         0, 0, 1];

    T_inv = [1/scale_factor, 0, centroid(1);
         0, 1/scale_factor, centroid(2);
         0, 0, 1];

end