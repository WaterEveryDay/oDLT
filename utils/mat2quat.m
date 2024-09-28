function q = mat2quat(A, scalar_first)
% Sebastien Henry
% transform SO(3) rotation matrix to quaternion
    theta = acos((trace(A) - 1) / 2);
    e = 1/(2 * sin(theta)) * [A(2,3) - A(3,2);
                              A(3,1) - A(1,3);
                              A(1,2) - A(2,1)];

    if scalar_first
        q = [cos(theta/2); e*sin(theta/2)];
    else
        q = [e*sin(theta/2); cos(theta/2)];
    end
end
