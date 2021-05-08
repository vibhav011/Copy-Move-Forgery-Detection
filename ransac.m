function [T, x0, inliers] = ransac(matches, N, k, beta)
    T = zeros(2, 2);
    A = zeros(k, 2);
    B = zeros(k, 2);
    x0 = zeros(2, 1);
    x1 = zeros(2, 1);
    x2 = zeros(2, 1);
    best_matched = 0;
    n = size(matches,1);
    
    for i = 1:N
        ids = randsample(n, k);
        
        A(:,:) = matches(ids, 1, :);
        B(:,:) = matches(ids, 2, :);
        muA = mean(A, 1);
        muB = mean(B, 1);
        
        A = A - muA;
        B = B - muB;
        
        cur_T = (B' * B) \ (B' * A);
        cur_x0 = muA' - cur_T * muB';
        
        num_inliers = 0;
        cur_inliers = zeros(n, 2, 2);
        for i2 = 1:n
            x1(:,:) = matches(i2, 1, :);
            x2(:,:) = matches(i2, 2, :);
            if norm(x1 - cur_T*x2 - cur_x0) < beta
                num_inliers = num_inliers + 1;
                cur_inliers(num_inliers, 1, :) = matches(i2, 1, :);
                cur_inliers(num_inliers, 2, :) = matches(i2, 2, :);
            end
        end
        
        if num_inliers > best_matched
            best_matched = num_inliers;
            inliers = cur_inliers(1:num_inliers, :, :);
            T = cur_T;
            x0 = cur_x0;
        end
    end
    
end