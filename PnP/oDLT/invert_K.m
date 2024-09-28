function K_inv = invert_K(K)
    dx = K(1, 1); dy = K(2, 2);
    cx = K(1, 3); cy = K(2, 3);
    a = K(1, 2);
    K_inv = [
        1/dx, -a /(dx*dy), (a*cy - dy*cx)/(dx*dy);
        0, 1/dy, -cy/dy;
        0, 0, 1];
end

