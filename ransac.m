function [T, x0, inliers] = ransac(matches, N, k, beta)
    T = zeros(2, 2);
    A = zeros(k, 2);
    B = zeros(k, 2);
    x0 = zeros(2, 1);
    x1 = zeros(2, 1);
    x2 = zeros(2, 1);
    best_matched = 0;
    n = size(matches,1);
    
    for i = 1:N % N=number of iterations for convergence of RANSAC - pre-decided 
        ids = randsample(n, k);
        
        A(:,:) = matches(ids, 1, :); % feature vectors of a set of k keypoints 
        B(:,:) = matches(ids, 2, :); % feature vectors of a the set of matches of the above k keypoints 
        muA = mean(A, 1); 
        muB = mean(B, 1);
        
        A = A - muA;
        B = B - muB;
        
        cur_T = (B' * B) \ (B' * A); % best estimate for T using the present set of k keypoints 
        cur_x0 = muA' - cur_T * muB'; % best estimate for x0 using the present set of k keypoints 
        
        num_inliers = 0; % counting number of inliers based on the above estimate 
        cur_inliers = zeros(n, 2, 2); 
        for i2 = 1:n
            x1(:,:) = matches(i2, 1, :); % feature vector of i2^{th} keypoint 
            x2(:,:) = matches(i2, 2, :); % feature vector of the i2^{th} keypoint's matching keypoint 
            if norm(x1 - cur_T*x2 - cur_x0) < beta % if cost function is less than the threshold beta
                num_inliers = num_inliers + 1; % number of inliers incremented 
                cur_inliers(num_inliers, 1, :) = matches(i2, 1, :); % recording the inliers and their matches
                cur_inliers(num_inliers, 2, :) = matches(i2, 2, :);
            end
        end
        
        if num_inliers > best_matched  % if present estimate gives more inliers than the best estimate till now, its results are recorded and values are updated 
            best_matched = num_inliers;
            inliers = cur_inliers(1:num_inliers, :, :);
            T = cur_T;
            x0 = cur_x0;
        end
    end
    
end
