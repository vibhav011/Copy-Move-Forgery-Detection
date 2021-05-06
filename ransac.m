function [T, x0, inliers] = ransac(matches, N, k, beta)
    T = zeros(2, 2);
    x0 = zeros(2, 1);
    best_matched = 0;
    n = size(matches);
    
    for i = 1:N
        ids = randsample(n, k);
        sum1 = zeros(2, 2);
        sum2 = zeros(2, 2);
        sum3 = zeros(2, 1);
        sum4 = zeros(2, 1);
        cur_T = T;
        cur_x0 = x0;
        error = 0;
        
        for j = 1:k
            sum1 = sum1 + matches(ids(j), 1) * matches(ids(j), 2)';
            sum2 = sum2 + matches(ids(j), 2) * matches(ids(j), 2)';
            sum3 = sum3 + matches(ids(j), 2)';
            sum4 = sum4 + matches(ids(j), 1);
            error = error + norm(matches(ids(j), 1) - cur_T*matches(ids(j), 2) - cur_x0)^2;
        end
        
        while error > 0.1
            cur_T = (sum1 - cur_x0*sum3) / sum2;
            cur_x0 = (sum4 - cur_T*sum3') ./ k;
            
            error = 0;
            for j = 1:k
                error = error + norm(matches(ids(j), 1) - cur_T*matches(ids(j), 2) - cur_x0)^2;
            end
        end
        
        num_inliers = 0;
        cur_inliers = [];
        for i2 = 1:n
            if norm(matches(i2, 1) - cur_T*matches(i2, 2) - cur_x0) < beta
                cur_inliers = [cur_inliers [matches(i2, 1) matches(i2, 2)]];
                num_inliers = num_inliers + 1;
            end
        end
        
        if num_inliers > best_matched
            best_matched = num_inliers;
            inliers = cur_inliers;
            T = cur_T;
            x0 = cur_x0;
        end
    end
    
end