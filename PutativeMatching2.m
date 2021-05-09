% X is 2 x n matrix of SIFT datapoints
% features is k x n matrix where n is the number of data points and k is dimension of each point
% pairs is a 2 x m matrix having m pairs of indices
function pairs = PutativeMatching2(X, features, alpha, offset)

[~, n] = size(features);
num_nbrs = min(n, (offset+1)^2+2);
num_out = round(alpha*n);

MdlKDT = KDTreeSearcher(features');
IdxKDT = knnsearch(MdlKDT, features', 'K', num_nbrs, 'Distance', 'euclidean');

dist = [];
pairs = [];

for i=1:n
    for j=1:num_nbrs
        l = IdxKDT(i, j);
        if abs(X(1, l)-X(1, i)) <= offset && abs(X(2, l)-X(2, i)) <= offset
            continue;
        end
        
        pairs = [pairs [i; l]];
        dist = [dist norm(features(:, l) - features(:, i))];
    end
end

[~, sortIdx] = sort(dist);
pairs = pairs(:, sortIdx);
pairs = pairs(:, 1:num_out);

end