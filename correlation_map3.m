function [map1, map2] = correlation_map3(img, W, transform, sz, x_min, x_max, y_min, y_max)
% outputs preliminary estimation for the patches as map1, map2 
[row,col]=size(img);
map1=zeros(row,col);
map2=zeros(row,col);

% x_min, x_max and y_min, y_max give the zoomed in window around the copied patches 
for x=x_min:x_max
    for y=y_min:y_max
        x_ = round(transform([x;y]));
        tx = x_(1);
        ty = x_(2);
        
        if tx <= sz || tx > row-sz || ty <= sz || ty > col-sz
            continue;
        end
        
        % path of size 2sz*2sz considered for findinf correlation co-efficient 
        p1 = img(tx-sz:tx+sz, ty-sz:ty+sz);
        p2 = W(tx-sz:tx+sz, ty-sz:ty+sz);
        
        R = corrcoef([reshape(p1, [], 1) reshape(p2, [], 1)]);
        corr = R(1,2);
        
        map1(x,y) = corr;
        map2(tx,ty) = max(map2(tx,ty), corr);
    end
end
end
