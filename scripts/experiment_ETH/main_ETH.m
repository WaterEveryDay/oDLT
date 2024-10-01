clear all
close all
clc

addpath ../../PnP/oDLT/;
addpath ../../PnP/EPnP/
addpath ../../PnP/RPnP/func;
addpath ../../PnP/CPnP/CPnP/
addpath ../../PnP/OPnP;
addpath ../../utils/

addpath utils/
addpath ../../utils/

profile off

scene_list = ["courtyard", "delivery_area", "electro", "facade", "kicker", "lounge", "meadow", "office", "playground", "relief", "terrace", "terrains"];

% change with your local folder to ETH3D stereo
% may download it at https://www.eth3d.net/
folder_to_eth3d = "/Users/sebastienhenry/Documents/datasets/eth3d/stereo/";

% output table containing results
filename = 'out/output_table.txt';

% methods to compare
method_names = ["DLT", "nDLT", "nDLT+GN", ...
                "EPNP", "EPnP+GN", "CPnP", ...
                "RPnP", "OPnP", ...
                "oDLT", "oDLT+LOST"];

funs = {@pnp_dlt, @pnp_dlt_normalized, @pnp_dlt_normalized_gn, ...
    @pnp_epnp, @pnp_epnp_gn, @pnp_cpnp, ...
    @pnp_rpnp, @pnp_opnp, ...
    @pnp_odlt, @pnp_odlt_lost};

linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--"];
markerstyles = ['p', 's', 'd', '*', 'x', '^', 'v', '>', '<', 'o'];
colors = colormap(turbo(10));

% remove DLT
method_names = method_names(2:end);
funs = funs(2:end);
linestyles = linestyles(2:end);
markerstyles = markerstyles(2:end);
colors = colors(2:end,:);

table_header(filename, method_names)
for scene_id = 1:length(scene_list)
    scene = scene_list(scene_id);
    type = "dslr"; % dslr: few views high res, rig: many views low res
    folder = folder_to_eth3d+scene+"/"+type+"_calibration_undistorted/";
    

    %% opening camera file
    camfile = fopen(folder+"/cameras.txt");

    % does not read the first 3 lines
    fgetl(camfile);
    fgetl(camfile);
    fgetl(camfile);

    % Define the format of the data
    formatSpec = '%d %s %d %d %f %f %f %f';

    % Read the data using textscan
    data = textscan(camfile, formatSpec);

    camparams = table(data{:}, 'VariableNames', {'CAMERA_ID', 'MODEL', 'WIDTH', 'HEIGHT', 'fx', 'fy', 'cx', 'cy'});
    nCameras = length(data{1});

    fclose(camfile);

    %% opening images file
    imagefile = fopen(folder+"/images.txt");
    fgetl(imagefile);
    fgetl(imagefile);
    fgetl(imagefile);
    formatSpec = '# Number of images: %d';
    nImages = fscanf(imagefile, formatSpec);

    % to finish the line
    fgetl(imagefile);

    CamData.f = zeros(1,nImages);
    CamData.k1 = zeros(1,nImages);
    CamData.k2 = zeros(1,nImages);
    CamData.T = zeros(3,3,nImages);
    CamData.r = zeros(3,nImages);
    camCalibMat = zeros(3, 3, nImages);

    % camera pose covariance
    SigmaRt = zeros(6, 6, nImages);

    % 2d point covariance on image plane (not pixel coordinates)
    Sigmaxx = zeros(2, 2, nImages);

    img_meas = cell(1, nImages);
    for i = 1:nImages
        formatSpec = '%d %f %f %f %f %f %f %f %d';
        line = fscanf(imagefile, formatSpec);
        image_ID = line(1);

        % quaternion
        Q = line(2:5);

        % translation
        T = line(6:8);

        camera_ID = line(9);

        fxi = camparams{camparams.CAMERA_ID==camera_ID,"fx"};
        fyi = camparams{camparams.CAMERA_ID==camera_ID,"fy"};
        cxi = camparams{camparams.CAMERA_ID==camera_ID,"cx"};
        cyi = camparams{camparams.CAMERA_ID==camera_ID,"cy"};

        camCalibMat(:,:,image_ID) = [fxi, 0, cxi; 0, fyi, cyi; 0, 0, 1];
        Sigmaxx(:,:,image_ID) = 1.5.^2*diag([1/fxi; 1/fyi].^2);

        CamData.T(:,:,image_ID) = quat2mat(Q, true)'; % R_CtoI
        CamData.r(:,image_ID) = -T; % TO VERIFY -posCam_C
        fgetl(imagefile);

        formatSpec2 = '%f %f %f';
        nextline =  fgetl(imagefile);
        elements = textscan(nextline, formatSpec2);
        img_meas(image_ID) = {[elements{:}]};
    end
    fclose(imagefile);

    %% opening point files (and tracks)
    pointfile = fopen(folder+"/points3D.txt");
    fgetl(pointfile);
    fgetl(pointfile);

    formatSpec = '# Number of points: %d';
    nPoints = fscanf(imagefile, formatSpec);
    % to finish the line
    fgetl(pointfile);

    Track.viewIds = [];
    Track.meas = [];
    Tracks = repmat(Track, nPoints, 1);

    % similarily to Tracks, construct a view structure that store the points
    % and the measurements
    % this will speed up the search for PnP (at a cost of more memory)
    View.pointIds = [];
    View.meas = [];
    Views = repmat(View, nImages, 1);

    % ETH3D dataset SfM solution (NOT scanner ground truth)
    x3dSol = zeros(3, nPoints);
    colorAll = zeros(3, nPoints);

    for ii = 1:nPoints
        formatSpec = '%d %f %f %f %d %d %d %f,';
        line = fscanf(pointfile, formatSpec);
        POINT3D_ID = line(1);

        % ETH3D dataset SfM solution (NOT scanner ground truth)
        X = line(2:4);
        x3dSol(:,ii) = X;

        % ETH3D dataset colors
        RGB = line(5:7);
        colorAll(:,ii) = RGB;

        ERROR = line(8);

        % IMAGE_ID POINT2D_IDX
        formatSpec2 = '%d %d';
        nextline =  fgetl(pointfile);
        elements = textscan(nextline, formatSpec2);
        [~, ia] = unique(elements{1}(:));
        elements(1) = {elements{1}(ia)};
        elements(2) = {elements{2}(ia)};

        for jj = 1:length(elements{1})
            VIEW_ID = elements{1}(jj);
            POINT2D_IDX = elements{2}(jj);
            K = camCalibMat(:,:,VIEW_ID);
            
            % transform pixel measurement to homogeneous coordinates
            meas = K\[img_meas{VIEW_ID}(POINT2D_IDX+1,1:2)';1];
            Tracks(ii).viewIds = [Tracks(ii).viewIds; VIEW_ID];
            Tracks(ii).meas = [Tracks(ii).meas; meas(1:2)'];

            Views(VIEW_ID).pointIds = [Views(VIEW_ID).pointIds; ii];
            Views(VIEW_ID).meas = [Views(VIEW_ID).meas; meas(1:2)'];
        end
    end
    fclose(pointfile);

    % create empty camera poses to be estimated
    camRots = zeros(3, 3, nImages); % R_ItoC
    camTrans = zeros(3, nImages); % - R_ItoC * posCam_I

    % this stores whether a camera is estimated
    isCamEstimated = boolean(zeros(nImages, 1));

    % this will store the number of 3D reconstructed points in each view
    nReconstructedPointsInViews = zeros(nImages, 1);

    % the estimate 3d points are set to not seen
    x3d = zeros(nPoints, 3);
    x3dCovs = zeros(3, 3, nPoints);
    is3dReconstructed = boolean(zeros(1, nPoints));
    nCamViewingPoints = zeros(1, nPoints);

    %% Construct a matrix counting how many common keypoints camera jj and kk see

    % a matrix that computes how many common points each camera see
    % this will come handy to decide what seed we want to choose
    % these are stored on the upper diagonal
    commonMap = zeros(nImages);
    matchedPoints = cell(nImages);

    for ii = 1:nPoints
        commonViews = Tracks(ii).viewIds;
        meas = Tracks(ii).meas;
        for jj = 1:length(commonViews)
            for kk = (jj+1):length(commonViews)
                v1 = jj;
                v2 = kk;
                if commonViews(v1) > commonViews(v2)
                    v1old = v1;
                    v1 = v2;
                    v2 = v1old;
                end
                idx1 = commonViews(v1);
                idx2 = commonViews(v2);
                commonMap(idx1,idx2) = commonMap(idx1,idx2) + 1;
                if commonMap(idx1, idx2) == 1
                    matchedPoints{idx1, idx2} = zeros(0,4);
                end
                matchedPoints{idx1, idx2} = [matchedPoints{idx1, idx2}; [meas(v1,:), meas(v2,:)]];
            end
        end
    end
    
    % profile on
    % compare_eth
    p = profile("info");
    % profile off
end