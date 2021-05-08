% X is 2 x n matrix of SIFT datapoints
% features is k x n matrix where n is the number of data points and k is dimension of each point
% pairs is a 2 x m matrix having m pairs of indices
function pairs = PutativeMatching(X, features, epsilon, offset)

[~, n] = size(features);
num_nbrs = min(n, (offset+1)^2+2);

MdlKDT = KDTreeSearcher(features');
IdxKDT = knnsearch(MdlKDT, features', 'K', num_nbrs, 'Distance', 'euclidean');

pairs = [];

for i=1:n
    cur_l = 0;
    for j=1:num_nbrs
        l = IdxKDT(i, j);
        if abs(X(1, l)-X(1, i)) <= offset && abs(X(2, l)-X(2, i)) <= offset
            continue;
        end
        
        if cur_l == 0
            cur_l = l;
        else
            if norm(features(:, cur_l) - features(:, i)) < epsilon*norm(features(:, l) - features(:, i))
                pairs = [pairs [i; cur_l]];
            end
            break;
        end
    end
end