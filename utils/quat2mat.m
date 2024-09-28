function A = quat2mat(q, scalar_first)
% Sebastien Henry
% transform quaternion to SO(3) rotation matrix
if scalar_first
    qs = q(1);
    qv = q(2:4);
else
    qs = q(4);
    qv = q(1:3);

end
qv = reshape(qv, 3, 1);
A = (qs^2 - qv'*qv) * eye(3) - 2 * qs * cross_mat(qv) + 2 * (qv * qv');
end
