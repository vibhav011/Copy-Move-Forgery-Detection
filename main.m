clear;
close all;

H=256;
W=256;
img_orig=imread('im9_t.bmp');
[H_orig, W_orig, ~] = size(img_orig);
img = imresize(img_orig,[H,W]);
img = im2double(rgb2gray(img));
s1 = W_orig/W;
s2 = H_orig/H;

[X,features] = SIFT(img, H, W);
imshow(img_orig);
hold on;
scatter(X(2, :)*s1, X(1, :)*s2, 'r.');
hold off;

pairs = PutativeMatching(X, features, 0.5, 5);
figure;
imshow(img_orig);
hold on;
plot(s1*[X(2,pairs(1,:)); X(2,pairs(2,:))], s2*[X(1,pairs(1,:)); X(1,pairs(2,:))], 'g');
hold off;

matches = zeros(size(pairs,2), 2, 2);
matches(:, 1, :) = X(:,pairs(1,:))';
matches(:, 2, :) = X(:,pairs(2,:))';
[T, x0, inliers] = ransac(matches, 500, 3, 3);
figure;
imshow(img_orig);
hold on;
plot(s1*[inliers(:, 1, 2)'; inliers(:, 2, 2)'], s2*[inliers(:, 1, 1)'; inliers(:, 2, 1)'], 'g');
hold off;

Tr = @(x) T*x+x0;
map = correlation_map(img,0.99, Tr);
new_img = zeros(H, W);
new_img(map == 1) = 128;
new_img(map == 2) = 255;
imshow(uint8(new_img));