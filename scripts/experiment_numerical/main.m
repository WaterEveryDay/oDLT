close all
clear all
clc
%%
addpath ../../PnP/oDLT/;
addpath ../../PnP/EPnP/
addpath ../../PnP/RPnP/func;
addpath ../../PnP/CPnP/CPnP/
addpath ../../PnP/OPnP;
addpath ../../utils/

rng(1)

%% parameters
% camera rotation
theta = 0*pi/180;
gamma = 0*pi/180;
R_I2C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)] ...
    * [cos(gamma), sin(gamma), 0; -sin(gamma), cos(gamma), 0; 0, 0, 1];

% camera position
r = [0; 0; 0];

% camera intrinsics
dx = 0.001250;
dy = dx;
K = [1/dx, 0, 320; 0, 1/dy, 240; 0, 0, 1]; %eye(3);

verbose = true;

% Monte Carlo
nmc = 100;

mode = "centered";
precise_timing = false;

% default parameters
n_meas = 50;
sig_u = 1;
outlier_percentage = 0;

% parameter changing, select
var_changing = "n";

n_meas_vec = 10:10:100;
n_meas_vec_timing = floor(logspace(1, 3.3, 20));
% n_meas_vec = [10, 100, 500, 1000];
sig_u_vec = [0.001, 1:5];
outlier_percentage_vec = [0, 5, 10, 15, 20, 25];

% calibrated experiment?
calibrated = true;

switch var_changing
    case "n"
        label = "number of measurements";
        n_changing = length(n_meas_vec);
        vec_changing = n_meas_vec;
    case "noise"
        label = "pixel noise";
        n_changing = length(sig_u_vec);
        vec_changing = sig_u_vec;
    case "outlier"
        label = "outlier percentage";
        n_changing = length(outlier_percentage_vec);
        vec_changing = outlier_percentage_vec;
    case "timing"
        label = "number of measurements";
        n_changing = length(n_meas_vec_timing);
        n_meas_vec = n_meas_vec_timing;
        vec_changing = n_meas_vec_timing;
        var_changing = "n";

    otherwise
        error("need to select a valid case")
end


%%
Sig_uu = diag([sig_u^2, sig_u^2, 0]);

%% define method names and functions to call
if calibrated
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
else
    idx = [2, 9];
    method_names = ["DLT", "nDLT (uncal)", "nDLT+GN", ...
    "EPNP", "EPnP+GN", "CPnP", ...
    "RPnP", "OPnP", ...
    "oDLT (uncal) (ours)", "oDLT+LOST (ours)"];

funs = {@pnp_dlt, @pnp_dlt_normalized_uncalibrated, @pnp_dlt_normalized_gn, ...
    @pnp_epnp, @pnp_epnp_gn, @pnp_cpnp, ...
    @pnp_rpnp, @pnp_opnp, ...
    @pnp_odlt_uncalibrated, @pnp_odlt_lost};

linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--"];
markerstyles = ['p', 's', 'd', '*', 'x', '^', 'v', '>', '<', 'o'];
colors = colormap(turbo(10));


method_names = method_names(idx);
funs = funs(idx);
linestyles = linestyles(idx);
markerstyles = markerstyles(idx);
colors = colors(idx, :);

n_methods = length(method_names);
end

%% define bounding boxes for spawning 3D points
switch mode
    case "centered"
        bbox = [-2, -2, 4; 2, 2, 8];
    case "uncentered"
        bbox = [1, 1, 4; 2, 2, 8];
    case "own1"
        %bbox = [-5, -5, 10; 5, 5, 15];
        bbox = [-5, -5, 10; 5, 5, 15];
end

%% define pixel boundaries
ubox = [0, 0; 640, 480];

%%
% a one number per Monte-Carlo to plot easily
aggregate_err_rot_all = zeros(n_changing, n_methods);
aggregate_err_rot_95_percentile_all = zeros(n_changing, n_methods);
aggregate_err_pos_all = zeros(n_changing, n_methods);
aggregate_err_reproj_all = zeros(n_changing, n_methods);
aggregate_err_fx_all = zeros(n_changing, n_methods);
aggregate_err_fy_all = zeros(n_changing, n_methods);
aggregate_err_cx_all = zeros(n_changing, n_methods);
aggregate_err_cy_all = zeros(n_changing, n_methods);
aggregate_time_all = zeros(n_changing, n_methods);

% The true projection from 3D to pixel
H = K*R_I2C*[eye(3), -r];
disp(R_I2C)

% iterate over the number of measurements.
for ii = 1:n_changing
    switch var_changing
        case "n"
            n_meas = n_meas_vec(ii);
        case "noise"
            sig_u = sig_u_vec(ii);
            Sig_uu = diag([sig_u^2, sig_u^2, 0]);
        case "outlier"
            outlier_percentage = outlier_percentage_vec(ii);
    end
    
    err_rot_all = zeros(nmc, n_methods);
    err_pos_all = zeros(3, nmc, n_methods);
    err_reproj_all = zeros(nmc, n_methods);
    err_fx_all = zeros(nmc, n_methods);
    err_fy_all = zeros(nmc, n_methods);
    err_cx_all = zeros(nmc, n_methods);
    err_cy_all = zeros(nmc, n_methods);
    time_all = zeros(nmc, n_methods);

    % iterate over the monte carlo
    for i = 1:nmc
        % generate 3D points randomly
        X = unifrnd(repmat(bbox(1,:)', 1, n_meas), ...
            repmat(bbox(2,:)',1, n_meas), 3, n_meas);

        % project 3D points to pixel
        U = H*[X; ones(1, n_meas)];
        U = U./U(3,:);

        % make pixel measurement noisy
        Utilde = U + mvnrnd([0; 0; 0], Sig_uu, n_meas)';

        % outlier ratio
        % corrupt first measurements according to outlier ratio
        n_outlier = floor(outlier_percentage / 100 * n_meas);
        Utilde(1:2, 1:n_outlier) = unifrnd(repmat(ubox(1,:)', 1, n_outlier), ...
            repmat(ubox(2,:)',1, n_outlier), 2, n_outlier);
        % compute pose for each method
        for method_id = 1:n_methods
            n_args_in = nargin(funs{method_id});
            n_args_out = nargout(funs{method_id});
            
            switch n_args_in
                case 2
                    args = {X', Utilde'};
                case 3
                    args = {X', Utilde', K};
               case 4
                    args = {X', Utilde', K};
            end
            
            methodtimer = tic;
            switch n_args_out
                case 2
                    [R_hat, t_hat] = funs{method_id}(args{:});
                    K_hat = K;
                case 3
                    [R_hat, t_hat, K_hat] = funs{method_id}(args{:});
            end
            runtime = toc(methodtimer);
            r_hat = -R_hat'*t_hat;

            if precise_timing
                f = @() funs{method_id}(args{:});
                runtime = timeit(f);
            end

            err_rot_all(i, method_id) = err_DCM(R_hat, R_I2C);
            err_pos_all(:, i, method_id) = r_hat - r;
            Hhat = K_hat*R_hat*[eye(3), -r_hat];
            err_reproj_all(i, method_id) = reprojection_error_using_matrix(Utilde, X, Hhat);
            time_all(i, method_id) = runtime;

            err_fx_all(i, method_id) = K(1,1) - K_hat(1,1);
            err_fy_all(i, method_id) = K(2,2) - K_hat(2,2);
            err_cx_all(i, method_id) = K(1,3) - K_hat(1,3);
            err_cy_all(i, method_id) = K(2,3) - K_hat(2,3);

        end
    end
    
    % aggregate results over all monte carlo samples for each method
    aggregate_err_rot_all(ii, :) = sqrt(mean(err_rot_all.^2, 1))*180/pi;
    aggregate_err_rot_95_percentile_all(ii,:) = prctile(err_rot_all, 95, 1)*180/pi;
    for method_id = 1:n_methods
        aggregate_err_pos_all(ii, method_id) = sqrt(mean(vecnorm(err_pos_all(:,:,method_id), 2).^2));
    end
    aggregate_err_reproj_all(ii,:) = mean(err_reproj_all, 1);
    aggregate_time_all(ii,:) = mean(time_all, 1);

    aggregate_err_fx_all(ii,:) = sqrt(mean(err_fx_all.^2, 1));
    aggregate_err_fy_all(ii,:) = sqrt(mean(err_fy_all.^2, 1));
    aggregate_err_cx_all(ii,:) = sqrt(mean(err_cx_all.^2, 1));
    aggregate_err_cy_all(ii,:) = sqrt(mean(err_cy_all.^2, 1));

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
        
        if ~calibrated
            disp("### fx RMSE")
            for method_id = 1:n_methods
                fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_fx_all(ii, method_id));
            end

            disp("### fy RMSE")
            for method_id = 1:n_methods
                fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_fy_all(ii, method_id));
            end

            disp("### cx RMSE")
            for method_id = 1:n_methods
                fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_cx_all(ii, method_id));
            end

            disp("### cy RMSE")
            for method_id = 1:n_methods
                fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_err_cy_all(ii, method_id));
            end
        end




        disp("### run times (s)")
        for method_id = 1:n_methods
            fprintf('%-24s: %8.5f\n', method_names(method_id), aggregate_time_all(ii, method_id));
        end
        fprintf("\n \n")
    end
end

%% figures
figure
fontsize_ax = 16;
fontsize_label = 20;
set(gcf,'Position',[100 100 1600 300])

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
switch var_changing
    case "n"
        switch mode
            case "uncentered"
                ylim([0, 5])
        end
    case "noise"
        switch mode
            case "uncentered"
                ylim([0,10])
        end
end

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
switch var_changing
    case "n"
        switch mode
            case "uncentered"
                ylim([0, 0.5])
        end
    case "noise"
        switch mode
            case "uncentered"
                ylim([0,1])
        end
end

set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("pos RMSE (unit)", fontsize=fontsize_label)
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
switch var_changing
    case "n"
        switch mode
            case "centered"
                ylim([1, 3])
            case "uncentered"
                ylim([1, 5])
        end
    case "noise"
        switch mode
            case "uncentered"
                ylim([0,30])
        end

end

set(gca,'FontSize',fontsize_ax)
xlabel(label, fontsize=fontsize_label)
ylabel("mean reprojection error (pixel)", fontsize=fontsize_label)
xlim([min(vec_changing), max(vec_changing)])
box on;
set(gca,'linewidth',2)
grid on;


% time
t4 = nexttile;
for method_id = 1:n_methods
    if method_id > 1
    hold on
    end
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
exportgraphics(t4, "figs/runtimes.pdf")


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
    newAxes.XLabel.String = '';
    newAxes.YLabel.String = '';

    % Adjust the position of the copied axes to fill the figure
    newAxes.Position = [0.1 0.05 0.9 0.9];

    set(newAxes,'FontSize',fontsize_ax)
    set(newAxes,'linewidth',2)
    
    % Save the figure
    switch i
        case 4
            filename = 'figs/'+var_changing+'_'+mode+'_rotation_plot.pdf';
        case 3
            filename = 'figs/'+var_changing+'_'+mode+'_position_plot.pdf';
        case 2
            filename = 'figs/'+var_changing+'_'+mode+'_reprojection_plot.pdf';
        case 1
            filename = 'figs/'+var_changing+'_'+mode+'_runtime_plot.pdf';
    end
    set(fig, 'PaperSize', [figWidth/100 figHeight/100]);
    print(fig, filename, '-dpdf', '-vector', '-bestfit')
    % print(fig, filename, '-dpng', '')

    % saveas(fig, filename)
    
    % Close the figure
    close(fig)
end

legend(method_names);


%% graph for the legends
figure
set(gcf, 'Position', [500, 500, 1000, 200]);
tiledlayout(1,1, Padding="compact", TileSpacing="tight")
nexttile
hold on
for method_id = 1:n_methods
    hold on
    semilogy(nan, nan, ...
        "LineStyle",linestyles(method_id), ...
        "Marker", markerstyles(method_id), ...
        LineWidth=2, Color=colors(method_id,:));
    hold off
end
lgd = legend(method_names, "FontSize",12, "Location", "north",'Orientation','horizontal');
%lgd.NumColumns=2;
lgd.Color = [0.95, 0.95, 0.95];
lgd.Box = 'on';
ax = gca;
ax.Visible = 'off';
exportgraphics(gcf, "figs/legends.pdf")


%% noisy measurement figure
figure
hold on
scatter(U(1,:), U(2,:));
scatter(Utilde(1,:), Utilde(2,:), '*');
legend("truth", "meas")

%% functions
function err = error_geometric(U, X, Phat)
% U = 3xnmeas
% X = 4xnmeas
% H = camera matrix
sums = 0;
n_meas = size(U, 2);
U = U./U(3,:);
Uhat = Phat*X;
Uhat = Uhat./Uhat(3,:);
for i = 1:n_meas
    sums = sums + norm(U(:,i) - Uhat(:,i))^2;
end
err = sum(sqrt((Uhat(1,:)-U(1,:)).^2+(Uhat(2,:)-U(2,:)).^2))/n_meas;
%err = sums/n_meas;
%err = norm(U - Uhat)^2;
end

function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)
%clear all; close all; load reprojection_error_usingRT;
n=size(Xw,1);

P=A*[R,T];
Xw_h=[Xw,ones(n,1)];
Urep_=(P*Xw_h')';

%project reference points into the image plane
Urep=zeros(n,2);
Urep(:,1)=Urep_(:,1)./Urep_(:,3);
Urep(:,2)=Urep_(:,2)./Urep_(:,3);

%reprojection error
err_=sqrt((U(:,1)-Urep(:,1)).^2+(U(:,2)-Urep(:,2)).^2);
err=sum(err_)/n;
end


