close all
clear all
clc
vesta_pt_cloud = pcread("/Users/sebastienhenry/Documents/datasets/Astrovision/vesta_3145k.ply");
vesta_crater_table = readtable("/Users/sebastienhenry/Documents/datasets/vesta_crater_database/database_700m.xlsx");

%%
a_craters = longlat2cartesian(vesta_crater_table{:,'LONGITUDE_CIRCLE'}'*pi/180, vesta_crater_table{:,'LATITUDE_CIRCLE'}'*pi/180);
crater_threshold = vesta_crater_table{:,"DIAM_CIRCLE"} < 100000 & vesta_crater_table{:,"DIAM_CIRCLE"} > 700;
theta = pi/180*-150;
R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
a_craters = a_craters(:,crater_threshold);
a_craters = R*a_craters;
[p_craters, altitudes] = getAltitudesFromPointCloud(vesta_pt_cloud.Location', a_craters);


figure
% Display point cloud
h = pcshow(vesta_pt_cloud);
hold on

% Camera parameters
camPos = get(gca, 'CameraPosition');
camTgt = get(gca, 'CameraTarget');
camDir = camTgt - camPos;
camUp = get(gca, 'CameraUpVector');
clear gca

% Extract point cloud data
vesta_points = vesta_pt_cloud.Location;  % Assuming this is Nx3 array

% Compute depth of each point (distance to camera)
depth = dot(vesta_points - camPos, repmat(camDir, size(vesta_points, 1), 1), 2);

% Set a threshold to filter out background points (adjustable)
depth_threshold = prctile(depth, 100);  % Only keep the closest 90% points
foreground_idx = depth <= depth_threshold;


% Compute depth of each point (distance to camera)
depth_crater = dot(p_craters' - camPos, repmat(camDir, size(p_craters', 1), 1), 2);
foreground_idx_crater = depth_crater <= depth_threshold;
%%
figure
normValues = vecnorm(vesta_points, 2, 2);
% Plot the filtered point cloud
pcshow(vesta_points(foreground_idx, :), normValues(foreground_idx), 'MarkerSize', 10);
hold on
colormap turbo
% Plot the craters in white
% scatter3(p_craters(1,foreground_idx_crater), p_craters(2,foreground_idx_crater), p_craters(3,foreground_idx_crater), 100, 'w', 'filled')
hold off
axis off
grid off
set(gcf, 'Color', 'w'); 

%% load image
addpath 'utils'
dataset_name = "astrovision";
scene_name = "2011205_rc3b";
scene_name = "2011218_opnav_023";
scene_name = "2011260_opnav_003";
scene_path = "/Users/sebastienhenry/Documents/datasets/Astrovision/"+scene_name+"/";

% read calibration file
load(scene_path+"params.mat")
Params = cell2mat(struct2cell(paramlist));

% creating structures to for camera intrinsics
fx = Params(1,1);
fy = Params(1,2);
cx = Params(1,3);
cy = Params(1,4);
K = [fx, 0, cx; 0, fy, cy; 0, 0, 1];

d = dir(scene_path+"images");
filenames = cellfun(@(x) fullfile(char(scene_path+"images/"), x), {d.name}, 'UniformOutput', false);
filenames = ""+filenames(3:end);
is_im_rgb = true;

ground_truth_poses = read_ground_truth_poses_astrovision(scene_path+"poses.mat");

addpath ../functions/
addpath ../EPnP;
addpath ../OPnP;
addpath ../MLPnP;
addpath ../functions;
addpath ../RPnP/func;
addpath ../CPnP/CPnP/
addpath ../OPnP;
method_names = ["DLT", "nDLT", "nDLT+GN", ...
    "EPNP", "EPnP+GN", "CPnP", ...
    "RPnP", "OPnP", ...
    "oDLT (ours)", "oDLT+LOST (ours)"];

funs = {@pnp_dlt, @pnp_dlt_normalized, @pnp_dlt_normalized_gn, ...
    @pnp_epnp, @pnp_epnp_gn, @pnp_cpnp, ...
    @pnp_rpnp, @pnp_opnp, ...
    @pnp_odlt, @pnp_odlt_lost};

linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--"];
markerstyles = ['p', 's', 'd', '*', 'x', '^', 'v', '>', '<', 'o'];
colors = colormap(turbo(10));
n_methods = length(method_names);

verbose = true;
n_changing = length(filenames);
aggregate_err_rot_all = zeros(n_changing, n_methods);
aggregate_err_rot_95_percentile_all = zeros(n_changing, n_methods);
aggregate_err_pos_all = zeros(n_changing, n_methods);
aggregate_err_reproj_all = zeros(n_changing, n_methods);
aggregate_time_all = zeros(n_changing, n_methods);
aggregate_converged_all = zeros(n_changing, n_methods);

nmc = 100;
vec_changing = 1:n_changing;
var_changing = "image number";
for ii = 1:length(vec_changing)

im_idx = vec_changing(ii);
im = imread(filenames(im_idx));


R = ground_truth_poses(im_idx).rotation;
r = ground_truth_poses(im_idx).position;

l = R*(p_craters - r);

% prune craters on the other side of the planet
fov_prune = (r/norm(r))'*(p_craters./altitudes) > cosd(90);
l = l(:, fov_prune);
p_craters_in_fov = p_craters(:, fov_prune);
x = l./l(3,:);
u_craters = K*x;

% prune craters outside fov
fov_prune = u_craters(1,:) > 0 & u_craters(1,:) < 1024 & u_craters(2,:) > 0 & u_craters(2,:) < 1024;
p_craters_in_fov = p_craters_in_fov(:, fov_prune);
u_craters = u_craters(:,fov_prune);

Sig_u = 1^2*[1, 0, 0; 0, 1, 0; 0, 0, 0];
n_meas = size(u_craters,2);


figure
imshow(im)
hold on
scatter(u_craters(1,:), u_craters(2,:), 'red', '+', LineWidth=2)
hold off

    err_rot_all = zeros(nmc, n_methods);
    err_pos_all = zeros(3, nmc, n_methods);
    err_reproj_all = zeros(nmc, n_methods);
    time_all = zeros(nmc, n_methods);
    converged_all = zeros(nmc, n_methods);

    for i = 1:nmc
        u_craters_tilde = u_craters + mvnrnd([0;0;0],Sig_u, n_meas)';
        for method_id = 1:n_methods
            methodtimer = tic;
            try
                [R_hat, t_hat] = funs{method_id}(p_craters_in_fov', u_craters_tilde', K);
                converged_all(i, method_id) = 1;
            catch
                R_hat = NaN(3);
                t_hat = NaN(3,1);
                converged_all(i, method_id) = 0;
            end
            runtime = toc(methodtimer);
            r_hat = -R_hat'*t_hat;

            err_rot_all(i, method_id) = err_DCM(R_hat, R);
            err_pos_all(:, i, method_id) = r_hat - r;
            Hhat = K*R_hat*[eye(3), -r_hat];
            err_reproj_all(i, method_id) = reprojection_error_using_matrix(u_craters_tilde, p_craters_in_fov, Hhat);
            % err_reproj_all(i, method_id) = reprojection_error_usingRT(X', Utilde', R_hat, -R_hat*r_hat, K);
            time_all(i, method_id) = runtime;
        end
    end
    % aggregate results over all monte carlo samples for each method
    aggregate_err_rot_all(ii, :) = sqrt(mean(err_rot_all.^2, 1, "omitnan"))*180/pi;
    aggregate_err_rot_95_percentile_all(ii,:) = prctile(err_rot_all, 95, 1)*180/pi;
    for method_id = 1:n_methods
        aggregate_err_pos_all(ii, method_id) = sqrt(mean(vecnorm(err_pos_all(:,:,method_id), 2).^2, "omitnan"));
    end
    aggregate_err_reproj_all(ii,:) = mean(err_reproj_all, 1, "omitnan");
    aggregate_time_all(ii,:) = mean(time_all, 1);
    aggregate_converged_all(ii,:) = mean(converged_all, 1);

    if verbose
        fprintf("n_meas = %i\n", n_meas)
        disp("### rotation errors RMSE (deg)")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_rot_all(ii, method_id));
        end
        disp("### rotation errors 95th percentile (deg)")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_rot_95_percentile_all(ii, method_id));
        end

        disp("### pos RMSE (units)")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_pos_all(ii, method_id));
        end
        disp("### mean geometric error")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_reproj_all(ii, method_id));
        end
        disp("### convergence rate")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_converged_all(ii, method_id));
        end
        disp("### run times (s)")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_time_all(ii, method_id));
        end
        fprintf("\n \n")
    end
end

%%
figure
fontsize_ax = 16;
fontsize_label = 20;
label = "image number";
set(gcf,'Position',[100 100 1600 400])
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% rotation error
t1 = nexttile;
for method_id = 1:n_methods
    hold on
    semilogy(vec_changing, aggregate_err_rot_all(:, method_id), ...
        "LineStyle",linestyles(method_id), ...
        "Marker", markerstyles(method_id), ...
        LineWidth=2, Color=colors(method_id,:))
    hold off
end
ylim([0, 3])


set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("rotation RMSE (deg)", fontsize=fontsize_label)
xlim([min(vec_changing), max(vec_changing)])
box on;
set(gca,'linewidth',2)
grid on;

% position error
t2 = nexttile;
for method_id = 1:n_methods
    hold on
    semilogy(vec_changing, aggregate_err_pos_all(:, method_id), ...
        "LineStyle",linestyles(method_id), ...
        "Marker", markerstyles(method_id), ...
        LineWidth=2, Color=colors(method_id,:))
    hold off
end
ylim([0, 50])

set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("pos RMSE (km)", fontsize=fontsize_label)
xlim([min(vec_changing), max(vec_changing)])
box on;
set(gca,'linewidth',2)
grid on;

% reproj error
t3 = nexttile;
for method_id = 1:n_methods
    hold on
    semilogy(vec_changing, aggregate_err_reproj_all(:, method_id), ...
        "LineStyle",linestyles(method_id), ...
        "Marker", markerstyles(method_id), ...
        LineWidth=2, Color=colors(method_id,:))
    hold off
end
hold on
set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("mean reprojection error (pixel)", fontsize=fontsize_label)
xlim([min(vec_changing), max(vec_changing)])
box on;
set(gca,'linewidth',2)
grid on;
ylim([1,3])


% time
t4 = nexttile;
for method_id = 1:n_methods
    loglog(vec_changing, aggregate_time_all(:, method_id)*1000, ...
        "LineStyle",linestyles(method_id), ...
        "Marker", markerstyles(method_id), ...
        LineWidth=2, Color=colors(method_id,:))
    hold off
end
% ylim([0, 12])
hold on
set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("mean runtime (ms)", fontsize=fontsize_label)
xlim([min(vec_changing), max(vec_changing)])
box on;
set(gca,'linewidth',2)
grid on;


%%

% Get the tile objects
tiles = t.Children;

figWidth = 400;
figHeight = 300;
% Loop through each tile
for i = 1:length(tiles)
    % Create a new figure with consistent size
    fig = figure('Visible', 'off', 'Position', [100, 100, figWidth, figHeight]);
    
    % Copy the content of the tile to the new figure
    newAxes = copyobj(tiles(i), fig);
    %newAxes.XLabel.String = '';
    %newAxes.YLabel.String = '';

    % Adjust the position of the copied axes to fill the figure
    newAxes.Position = [0.1 0.1 0.9 0.9];

    set(newAxes,'FontSize',fontsize_ax)
    set(newAxes,'linewidth',2)
    
    % Save the figure
    switch i
        case 4
            filename = 'figs/'+var_changing+'_rotation_plot.pdf';
        case 3
            filename = 'figs/'+var_changing+'_position_plot.pdf';
        case 2
            filename = 'figs/'+var_changing+'_reprojection_plot.pdf';
        case 1
            filename = 'figs/'+var_changing+'_runtime_plot.pdf';
    end
    set(fig, 'PaperSize', [figWidth/80 figHeight/80]);
    print(fig, filename, '-dpdf', '-vector', '-bestfit')
    % print(fig, filename, '-dpng', '')

    % saveas(fig, filename)
    
    % Close the figure
    close(fig)
end

%%
figure
normValues = vecnorm(vesta_points, 2, 2);
% Plot the filtered point cloud
pcshow(vesta_points(foreground_idx, :), normValues(foreground_idx), 'MarkerSize', 10);
hold on
colormap turbo
% Plot the craters in white
% scatter3(p_craters(1,foreground_idx_crater), p_craters(2,foreground_idx_crater), p_craters(3,foreground_idx_crater), 100, 'w', 'filled')
hold off
axis off
grid off
set(gcf, 'Color', 'w'); 


highlight_view = 0;
for ii = 1:4:length(vec_changing)
    view_id = vec_changing(ii);
    hold on
    pose = rigidtform3d(ground_truth_poses(view_id).rotation', ground_truth_poses(view_id).position)  ;

    cam_color = "black";
    if view_id == highlight_view
        cam_color = "blue";
    end
    cam = plotCamera(AbsolutePose=pose,Opacity=0, Size=10, Color=cam_color);
    hold off
end

% ax2 = axes('Position',[0.65 0.65 0.28 0.28]);
% hold on
% scatter3(p_craters_in_fov(1,:), p_craters_in_fov(2,:), p_craters_in_fov(3,:), 100, 'w', 'filled')

% exportgraphics(gca, scene_list(scene_id)+"_reconstruction.pdf")
% % p = profile("info");
% % profile off
% disp("")


%%

function u = longlat2cartesian(long, lat)
u = [cos(long) .* cos(lat);
    sin(long) .* cos(lat);
    sin(lat)];
end

function [p_craters, altitudes] = getAltitudesFromPointCloud(vesta_locations, a_craters)
% Convert longitudes and latitudes to unit vectors
vesta_norms = vecnorm(vesta_locations, 2, 1);
vesta_unit = vesta_locations./vesta_norms;
p_craters = a_craters;

% Initialize altitudes vector
n = size(p_craters, 2);
altitudes = zeros(1,n);

% Loop through each unit vector to find the corresponding altitude
for i = 1:n
    % Find the closest point in the point cloud
    center_direction = a_craters(:, i);

    cos_angles = center_direction'*vesta_unit;

    % Find the index of the nearest point
    [~, nearest_idx] = max(cos_angles);

    % Get the altitude corresponding to the nearest point
    altitudes(i) = vesta_norms(nearest_idx); % Assuming altitude is in the 3rd row
    p_craters(:, i) = altitudes(i) * a_craters(:, i);
end
end


