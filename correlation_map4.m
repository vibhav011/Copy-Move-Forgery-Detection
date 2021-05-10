function [best_map1, best_map2] = correlation_map4(img, theta, alpha, p, q, sz, th)
[row,col]=size(img);

best_map1 = zeros(row, col);
best_map2 = zeros(row, col);
count = 0;
h = 1/3*ones(3,1);
H = h*h';

for dt = 0:0
    dt
    d2r = deg2rad(dt);
    T = alpha * [cos(theta+d2r) -sin(theta+d2r); sin(theta+d2r) cos(theta+d2r)]';
    x0 = p - T*q;
    
    transform = @(x) T*x+x0;
    tform = affine2d([T(2,2) T(1,2) 0;
                      T(2,1) T(1,1) 0;
                      x0(2) x0(1) 1]);
    
    W = imwarp(img, tform, 'OutputView',imref2d(size(img)));
    
    map1=zeros(row,col);
    map2=zeros(row,col);

    for x=1:row
        for y=1:col
            x_ = round(transform([x;y]));
            tx = x_(1);
            ty = x_(2);

            if tx <= sz || tx > row-sz || ty <= sz || ty > col-sz
                continue;
            end

            p1 = img(tx-sz:tx+sz, ty-sz:ty+sz);
            p2 = W(tx-sz:tx+sz, ty-sz:ty+sz);
%             p1 = filter2(H,p1);
%             p2 = filter2(H,p2);

            R = corrcoef([reshape(p1, [], 1) reshape(p2, [], 1)]);
            corr = R(1,2);

            map1(x,y) = corr;
            map2(tx,ty) = max(map2(tx,ty), corr);
        end
    end
    cur_count = sum(map1>th, 'all')
    if cur_count > count
        count = cur_count;
        best_map1 = map1;
        best_map2 = map2;
    end
end
end
