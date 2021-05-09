clear;
close all;

H=256;
W=256;
img_orig=imread('out_r-20_im8.bmp');
[H_orig, W_orig, ~] = size(img_orig);
% img = imresize(img_orig,[H,W]);
% img2 = imresize(img, 2*[H, W]);
img = im2single(rgb2gray(img_orig));

s1 = W_orig/W;
s2 = H_orig/H;
s1 = 1;
s2 = 1;

[X, features] = vl_sift(img);
X = round(X(1:2, :));
X = [X(2,:); X(1,:)];
features = double(features);
% [X,features] = SIFT_new(img, H, W);
figure;
imshow(img_orig);
drawnow;
hold on;
scatter(X(2, :)*s1, X(1, :)*s2, 'r.');
hold off;

pairs = PutativeMatching(X, features, 0.5, 5);
figure;
imshow(img_orig);
drawnow;
hold on;
plot(s1*[X(2,pairs(1,:)); X(2,pairs(2,:))], s2*[X(1,pairs(1,:)); X(1,pairs(2,:))], 'g');
hold off;

matches = zeros(size(pairs,2), 2, 2);
matches(:, 1, :) = X(:,pairs(1,:))';
matches(:, 2, :) = X(:,pairs(2,:))';
[T, x0, inliers] = ransac(matches, 100, 3, 6);
figure;
imshow(img_orig);
drawnow;
hold on;
plot(s1*[inliers(:, 1, 2)'; inliers(:, 2, 2)'], s2*[inliers(:, 1, 1)'; inliers(:, 2, 1)'], 'g');
hold off;
 
Tr = @(x) T*x+x0;
tform = affine2d([T(:,1)' 0;
                  T(:,2)' 0;
                  x0' 1]);
W = imwarp(img, tform);
figure;
imshow(img);
figure;
imshow(W);
% map = correlation_map(img,0.999, Tr);
% map = correlation_map2(img,0.7, Tr, 3);
% new_img = zeros(H_orig, W_orig);
% new_img(map == 1) = 255;
% new_img(map == 2) = 128;
% figure;
% imshow(uint8(new_img));