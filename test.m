%should write a path to the vl_setup.m file of vlfeat package
run('.\vlfeat-0.9.21\toolbox\vl_setup.m');


%prepare images

impath1 = fullfile('.\PA2_dataset\dataset_twoview\sfm01.jpg');
impath2 = fullfile('.\PA2_dataset\dataset_twoview\sfm02.jpg');

im1_rgb = imread(impath1);
im2_rgb = imread(impath2);

im1 = im2single(rgb2gray(im1_rgb));
im2 = im2single(rgb2gray(im2_rgb));

%calibration matrix
K = [ 1506.070 0 1021.043; 0 1512.965 698.031; 0 0 1 ];

%extract sift features
[f1, d1] = vl_sift(im1);
[f2, d2] = vl_sift(im2);

%calculate matching
[matches, scores] = vl_ubcmatch(d1, d2);

%slicing frames according to the matching
f1 = [f1(1:2, matches(1,:)); ones(1, size(matches, 2))];
f2 = [f2(1:2, matches(2,:)); ones(1, size(matches, 2))];


%%run RANSAC 5-points algorithm to eliminate outlier
iterations = 3000; % # of times to iterate the algorithm
t = 1;

finalE = []; % final Essential matrix calculated by RANSAC
finalInliers = [];
inlier_ratio = 0;


while t<iterations
    %ramdomly select five points
    idxs = randperm(size(matches,2));
    idxs = idxs(1:5);

    %slicing keypoints according to idxs
    %normalizing the image coordinates by calibration matrix
    f1_n = K\f1(:, idxs);
    f2_n = K\f2(:, idxs);

    %calculate essential matrices
    Es = calibrated_fivepoint(f1_n, f2_n);
  

    %reshape Es(9xN) into E(3x3xN) matrices
    E = [];
    for i=1:size(Es,2)
        E(:,:,i) = reshape(Es(:,i),[3,3]);
        
    end

    %for each essential matrix, calculate fundamental matrix and count the
    %number of inliers by calculating distance
    %so fianlly determine best Essential matrix which has the largest number of inliers: "finalE"
    for n=1:size(Es, 2)
        F = inv(K')*E(:,:,n)*inv(K);

        %calculating distance for each match keypoint
        for j=1:size(f1,2)
            d(j) = distance(F, f1(:,j), f2(:,j));
        end

        %count the number of points 
        inliers = find(d < 0.05); % threshhold 
        inliers_count = size(inliers, 2);
        total_count = size(f1, 2);

        if (inliers_count > size(finalInliers,2))
            finalInliers = inliers;
            finalE = E(:,:,n);
            inlier_ratio = inliers_count / total_count;
            disp(inlier_ratio);
        end


    end

    t = t+1;
end

[U, S, V] = svd(finalE); %decomposition of Essential matrix



W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 0 0];
u3 = U * Z * U';

%original camera matrix
P = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P = K * P;

%four possible relative camera matrix
P1 = K * horzcat(U*W*V', u3);
P2 = K * horzcat(U*W*V', -u3);
P3 = K * horzcat(U*W'*V', u3);
P4 = K * horzcat(U*W'*V', -u3);


worldX1 = GetWorldCoordinates(P, P1, f1, f2);
worldX2 = GetWorldCoordinates(P, P2, f1, f2);
worldX3 = GetWorldCoordinates(P, P3, f1, f2);
worldX4 = GetWorldCoordinates(P, P4, f1, f2);

color = [];


for i=1:size(f1,2)
    color_info = im1_rgb(round(f1(2,i)),round(f1(1,i)),:);
    color_info = reshape(color_info,[3,1]);
    color = horzcat(color,color_info);
end

SavePLY('X1.ply',[worldX1;color]);
SavePLY('X2.ply',[worldX2;color]);
SavePLY('X3.ply',[worldX3;color]);
SavePLY('X4.ply',[worldX4;color]);













