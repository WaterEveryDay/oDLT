function ground_truth_poses = read_ground_truth_poses_astrovision(ground_truth_path)
    % create a structure to store the ground truth
    load(ground_truth_path, 'translation', 'quaternion');
    ground_truth_poses(size(quaternion, 1)) =...
        struct( 'rotation', [], ...
                'position', [], ...
                'timestamp',[]);

    for i = 1:size(quaternion, 1)
        ground_truth_poses(i).timestamp = NaN;
        % convention on this dataset: R_C_to_I
        % such that x_C = R_I_to_C * x_I
        ground_truth_poses(i).rotation = quat2mat(quaternion(i,:)', true)';
        % convention on this dataset: t = -R_I_to_C*r_I
        % so we bring it back to r_I (inertial position of camera)
        ground_truth_poses(i).position = -ground_truth_poses(i).rotation'*translation(i,:)';
    end
end
    