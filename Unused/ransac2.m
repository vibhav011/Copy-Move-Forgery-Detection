function [T, alpha ,x0, inliers] = ransac2(matches, N, k, beta)
    [T1,x0,~] = ransac(matches,N,k,beta);
%    cur_T = eye(2,2);
%    T = eye(2, 2);
%     A = zeros(2, 1);
%     B = zeros(2, 1);
%     x = zeros(2, 1);
%     xm = zeros(2, 1);
    best_matched = 0;
    n = size(matches,1);
    ex=(T1*[1;0]);
    ey=(T1*[0;1]);
    alpha=(norm(ex)+norm(ey))/2;
    for i = 1:N % N=number of iterations for convergence of RANSAC - pre-decided 
        id = randsample(n, 1);
        x = matches(id, 1, :); % coordinates of a set of k keypoints 
        xm = matches(id, 2, :); % coordinates of matches of the above k keypoints 
        A = alpha*[x(1,1,1);x(1,1,2)];
        B = [xm(1,1,1)-x0(1,1);xm(1,1,2)-x0(2,1)];
%         cur_T = (B' * B) \ (B' * A); % best estimate for T using the present set of k keypoints 
%         cur_x0 = muA' - cur_T * muB'; % best estimate for x0 using the present set of k keypoints
        [U,~,V] = svd(A*B');
        cur_T = V*U';
        if(det(cur_T)<0)
            continue;
        end
         
        num_inliers = 0; % counting number of inliers based on the above estimate 
        cur_inliers = zeros(n, 2, 2); 
        for i2 = 1:n
            x1(:,:) = matches(i2, 1, :); % feature vector of i2^{th} keypoint 
            x2(:,:) = matches(i2, 2, :); % feature vector of the i2^{th} keypoint's matching keypoint 
            if norm(x1 - cur_T*x2 - x0) < beta || norm(x2 - cur_T*x1 - x0) < beta % if cost function is less than the threshold beta
                num_inliers = num_inliers + 1; % number of inliers incremented 
                cur_inliers(num_inliers, 1, :) = matches(i2, 1, :); % recording the inliers and their matches
                cur_inliers(num_inliers, 2, :) = matches(i2, 2, :);
            end
        end
        
        if num_inliers >= best_matched  % if present estimate gives more inliers than the best estimate till now, its results are recorded and values are updated 
            best_matched = num_inliers;
            inliers = cur_inliers(1:num_inliers, :, :);
            T = cur_T;
%             x0 = cur_x0;
        end
    end
    
end
