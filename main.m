clear;
close all;
rng(2);

H=256;
W=256;
img_orig=imread('images/out_s1.2_im9.bmp');
[H_orig, W_orig, ~] = size(img_orig);
% img = imresize(img_orig,[H,W]);
% img2 = imresize(img, 2*[H, W]);
img = im2single(rgb2gray(img_orig));

[X, features] = vl_sift(img);
X = round(X(1:2, :));
X = [X(2,:); X(1,:)];
features = double(features);
% [X,features] = SIFT_new(img, H, W);
figure;
imshow(img_orig);
drawnow;
hold on;
scatter(X(2, :), X(1, :), 'r.');
hold off;

pairs = PutativeMatching(X, features, 0.5, 5);
figure;
imshow(img_orig);
drawnow;
hold on;
plot([X(2,pairs(1,:)); X(2,pairs(2,:))], [X(1,pairs(1,:)); X(1,pairs(2,:))], 'g');
hold off;

matches = zeros(size(pairs,2), 2, 2);
matches(:, 1, :) = X(:,pairs(1,:))';
matches(:, 2, :) = X(:,pairs(2,:))';
[T1, ~, inliers] = ransac(matches, 200, 3, 3);
figure;
imshow(img_orig);
drawnow;
hold on;
plot([inliers(:, 1, 2)'; inliers(:, 2, 2)'], [inliers(:, 1, 1)'; inliers(:, 2, 1)'], 'g');
hold off;

[T, x0] = ransac3(img, inliers, T1);


% [U,S,V] = svd(T);
% S = mean([S(1,1) S(2,2)])*eye(2);
% T = U*S*V';

% T = [0.9505 -0.3409; 0.3596 0.9362];

Tr = @(x) T*x+x0;
tform = affine2d([T(2,2) T(1,2) 0;
                  T(2,1) T(1,1) 0;
                  x0(2) x0(1) 1]);


W = imwarp(img, tform, 'OutputView',imref2d(size(img)));
% figure;
% imshow(img);
figure;
imshow(W);
% map = correlation_map(img,0.999, Tr);
% map = correlation_map2(img,0.7, Tr, 3);
[map1, map2] = correlation_map3(img,W, Tr, 2);
mask1 = location(map1, 0.9, 0.9);
mask2 = location(map2, 0.9, 0.9);
mask = min(mask1+mask2, 1);

figure;
imshow(mask);