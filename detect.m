function out_mask = detect(img_orig)

H=256;
W=256;
% img_orig=imread('images/out_r-4_im15.bmp');
[H_orig, W_orig, ~] = size(img_orig);
% img = imresize(img_orig,[H,W]);
% img2 = imresize(img, 2*[H, W]);
img = im2single(rgb2gray(img_orig));

[X, features] = vl_sift(img);
X = round(X(1:2, :));
X = [X(2,:); X(1,:)];
features = double(features);
% [X,features] = SIFT_new(img, H, W);
% figure;
% imshow(img_orig);
% drawnow;
% hold on;
% scatter(X(2, :), X(1, :), 'r.');
% hold off;

pairs = PutativeMatching(X, features, 0.5, 5);
% figure;
% imshow(img_orig);
% drawnow;
% hold on;
% plot([X(2,pairs(1,:)); X(2,pairs(2,:))], [X(1,pairs(1,:)); X(1,pairs(2,:))], 'g');
% hold off;

matches = zeros(size(pairs,2), 2, 2);
matches(:, 1, :) = X(:,pairs(1,:))';
matches(:, 2, :) = X(:,pairs(2,:))';
[T1, ~, inliers] = ransac(matches, 200, 3, 3);

% figure;
% imshow(img_orig);
% drawnow;
% hold on;
% % plot([inliers(:, 1, 2)'; inliers(:, 2, 2)'], [inliers(:, 1, 1)'; inliers(:, 2, 1)'], 'g');
% scatter(inliers(:,1,2), inliers(:,1,1), 'r.');
% hold on;
% scatter(inliers(:,2,2), inliers(:,2,1), 'b.');
% hold off;

[T, x0] = ransac3(img, inliers, T1);

Tr = @(x) T*x+x0;
tform = affine2d([T(2,2) T(1,2) 0;
                  T(2,1) T(1,1) 0;
                  x0(2) x0(1) 1]);
               
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
% 
% 
W = imwarp(img, tform, 'OutputView',imref2d(size(img)));
% % figure;
% % imshow(img);
% figure;
% imshow(W);
% % map = correlation_map(img,0.999, Tr);
% map = correlation_map2(img,0.7, Tr, 3);
% [map1, map2] = correlation_map3(img,W, Tr, 2, 1, H_orig, 1, W_orig);
[map1, map2] = correlation_map3(img,W, Tr, 2, x_min, x_max, y_min, y_max);
% [map1, map2] = correlation_map4(img, theta, alpha, p, q, 2, 0.9);

[mask1, idx1] = location(map1, 0.45);
[mask2, idx2] = location(map2, 0.45);
out_mask = min(mask1+mask2, 1);

% figure;
% imshow(out_mask);
% drawnow;
end