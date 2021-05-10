% main driver code
clear;
close all;
rng(4); % seed 

%% pre-processing and reading 
H=256;
W=256;
img_orig=imread('images/im3_t.bmp'); % can replace with image of choice 
[H_orig, W_orig, ~] = size(img_orig);
img = im2single(rgb2gray(img_orig));

%% SIFT
[X, features] = vl_sift(img);
X = round(X(1:2, :));
X = [X(2,:); X(1,:)];
features = double(features);

% plotting the SIFT key-points
figure;
imshow(img_orig);
drawnow;
hold on;
scatter(X(2, :), X(1, :), 'r.');
hold off;

%% PUTATIVE MATCHING - threshold for (lowest distance)/(second lowest diatance)  = 0.7, offset = 5
pairs = PutativeMatching(X, features, 0.7, 5);
figure;
imshow(img_orig);
drawnow;
hold on;
plot([X(2,pairs(1,:)); X(2,pairs(2,:))], [X(1,pairs(1,:)); X(1,pairs(2,:))], 'g');
hold off;

%% RANSAC 
matches = zeros(size(pairs,2), 2, 2); % matches(i,j,k) contains the k^{th} (in {1,2}) co-ordinate of th  j^{th} (in {1,2}) component of the i^{th} pair of matches
matches(:, 1, :) = X(:,pairs(1,:))';
matches(:, 2, :) = X(:,pairs(2,:))';
[T1, ~, inliers] = ransac(matches, 200, 3, 3); % number of runs = 200, number of pairs in one iteration = 3, threshold on cost function = 3

figure;
imshow(img_orig);
drawnow;
hold on;
plot([inliers(:, 1, 2)'; inliers(:, 2, 2)'], [inliers(:, 1, 1)'; inliers(:, 2, 1)'], 'g');
hold off;

[T, x0] = ransac3(img, inliers, T1); % T gives the rotation matrix 

%% AFFINE TRANSFORMATION
Tr = @(x) T*x+x0; % affine transform that absorbs scaling, rotation and translation 
tform = affine2d([T(2,2) T(1,2) 0;
                  T(2,1) T(1,1) 0;
                  x0(2) x0(1) 1]);

W = imwarp(img, tform, 'OutputView',imref2d(size(img))); % affine transformation applied to entire image

figure;
imshow(W);

%% LOCATION AND CORRELATION
% zooming in into the neighbourhood of the matched patches (to generate correlation map) to increase speed and minimise number of  outliers 

in_ind = 1;
[num_in, ~, ~] = size(inliers);
a = zeros(num_in,2);
b = zeros(num_in,2);
a(:,:) = inliers(:,1,:);
b(:,:) = inliers(:,2,:);

if sum(vecnorm(b'-Tr(a'))) > sum(vecnorm(a'-Tr(b')))
    in_ind = 2;
end

x_min = min(inliers(:,in_ind,1));
x_max = max(inliers(:,in_ind,1));
y_min = min(inliers(:,in_ind,2));
y_max = max(inliers(:,in_ind,2));
diff_x = x_max-x_min;
diff_y = y_max-y_min;

x_min = max(1, x_min-2*diff_x);
x_max = min(H_orig, x_max+2*diff_x);
y_min = max(1, y_min-2*diff_y);
y_max = min(W_orig, y_max+2*diff_y);

% basic correlation check (produces noisy mask)
[map1, map2] = correlation_map3(img,W, Tr, 2, x_min, x_max, y_min, y_max);

% locating patches and smoothening the mask for both components of the pair of copied patches 
[mask1, idx1] = location(map1, 0.9); % correlation threshold = 0.9
[mask2, idx2] = location(map2, 0.9); % correlation threshold = 0.9
mask = min(mask1+mask2, 1); % both patch masks combined in one mask

figure;
imshow(mask);
