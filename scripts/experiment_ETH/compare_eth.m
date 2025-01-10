% profile on
% tested methods

n_methods = length(method_names);

% compute the average number of points per view
lengths = arrayfun(@(view) length(view.pointIds), Views); % Get lengths
n_avg = mean(lengths); % Compute the average

% iterate for every view
n_views = length(Views);
err_rot_all = zeros(n_views, n_methods);
err_pos_all = zeros(3, n_views, n_methods);
err_reproj_all = zeros(n_views, n_methods);
time_all = zeros(n_views, n_methods);
pose_all = zeros(3,4, n_views, n_methods);

% compute pose for each method
for method_id = 1:n_methods
    for view_id = 1:length(Views)
        point_ids = Views(view_id).pointIds;
        p_all = x3dSol(:, point_ids);
        u_all = Views(view_id).meas';
        Kreal = camCalibMat(:,:,view_id);
        K = Kreal; %eye(3);
        u_all = K*[u_all; ones(1, size(u_all, 2))];


        Sig_uu = diag([1, 1, 0]);
        R_W2C_sol = CamData.T(:,:,view_id);
        r_sol = R_W2C_sol'*CamData.r(:,view_id);

        
        nArgs = nargin(funs{method_id});
        rng(1)
        
        switch nArgs
            case 2
                args = {p_all', u_all'};
            case 3
                args = {p_all', u_all', K};
            case 4
                args = {p_all', u_all', K};
            case 5
                args = {p_all', u_all', K};
        end
        methodtimer = tic;
        [R_hat, t_hat] = funs{method_id}(args{:});
        runtime = toc(methodtimer);

        r_hat = -R_hat'*t_hat;
        err_rot_all(view_id, method_id) = err_DCM(R_hat, R_W2C_sol);
        err_pos_all(:, view_id, method_id) = r_hat - r_sol;
        Hhat = K*R_hat*[eye(3), -r_hat];
        err_reproj_all(view_id, method_id) = reprojection_error_using_matrix(u_all, p_all, Hhat);
        pose_all(:,:,view_id, method_id) = [R_hat, -R_hat*r_hat];
        time_all(view_id, method_id) = runtime;
    end
end
% a one number to plot easily
aggregate_err_rot_all = zeros(1, n_methods);
aggregate_err_pos_all = zeros(1, n_methods);
aggregate_err_reproj_all = zeros(1, n_methods);
aggregate_time_all = zeros(1, n_methods);

% aggregate results over all monte carlo samples for each method
aggregate_err_rot_all(1, :) = sqrt(mean(err_rot_all.^2, 1))*180/pi;
for method_id = 1:n_methods
    aggregate_err_pos_all(1, method_id) = sqrt(mean(vecnorm(err_pos_all(:,:,method_id), 2).^2));
end
aggregate_err_reproj_all(1,:) = mean(err_reproj_all, 1);
aggregate_time_all(1,:) = median(time_all, 1);

%%
disp("scene = " + scene)
disp("### rotation RMSE (deg)")
for method_id = 1:n_methods
    fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_rot_all(1, method_id));
end

disp("### pos RMSE (units)")
for method_id = 1:n_methods
    fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_pos_all(1, method_id));
end
disp("### mean geometric error (pixel)")
for method_id = 1:n_methods
    fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_reproj_all(1, method_id));
end
disp("### mean run times (ms)")
for method_id = 1:n_methods
    fprintf('%-24s: %8.5f\n', method_names(method_id), 1000*aggregate_time_all(1, method_id));
end
fprintf("\n \n")

write_latex_table(scene, method_names, aggregate_err_rot_all, aggregate_err_pos_all, aggregate_err_reproj_all, aggregate_time_all, n_avg, filename);

%% scatter reproj - time
% scatter(aggregate_time_all*1000, aggregate_err_reproj_all)
% xlabel("runtime (ms)", fontsize = 20)
% ylabel("reprojection error (pixel)", fontsize = 20)

% Your data (assuming aggregate_time_all and aggregate_err_reproj_all are defined)

aggregate_err_reproj_all = aggregate_err_reproj_all;

% Markers and labels (adjust based on the number of points)

figure;
hold on; % To plot all points on the same plot

% Plot each point individually with different markers and labels
for i = 1:length(aggregate_time_all)
    scatter(aggregate_time_all(i)*1000, aggregate_err_reproj_all(i), 100, ...
        'Marker', markerstyles(i), ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', method_names(i));
    %text(aggregate_time_all(i)*1000 + 0.1, aggregate_err_reproj_all(i), ['  ', method_names(i)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12); % Label each point
end
set(gca,'FontSize',14);
xlabel("runtime (ms)", fontsize = 20)
ylabel("mean reprojection error (pixel)", fontsize = 20)
[out1, out2] = legend(fontsize = 20);
for i = 1:n_methods
    out2(n_methods+i).Children.MarkerSize = 12;
    out2(n_methods+i).Children.LineWidth = 2;
end
lgd.ItemTokenSize = [60, 60];
xlim([0,15])
box on;
set(gca,'linewidth',2)
grid on;

exportgraphics(gca, "pareto.pdf")




%% plot camera
figure
highlight_view = n_views-3;
for view_id = 1:n_views
    hold on
    pose = rigidtform3d(pose_all(1:3, 1:3, view_id, end)',-pose_all(1:3, 1:3, view_id, end)'*pose_all(1:3, 4, view_id, end))  ;

    cam_color = "red";
    if view_id == highlight_view
        cam_color = "blue";
    end
    cam = plotCamera(AbsolutePose=pose,Opacity=0, Size=0.1, Color=cam_color);
    hold off
end

hold on
pcshow(x3dSol', colorAll'/256, BackgroundColor='w', MarkerSize=10)
hold off
xlim([-50, 50])
ylim([-50, 50])
exportgraphics(gca, scene_list(scene_id)+"_reconstruction.pdf")
% p = profile("info");
% profile off
disp("")
