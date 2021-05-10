function [map1, map2] = correlation_map3(img, W, transform, sz)
[row,col]=size(img);
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
        
        R = corrcoef([reshape(p1, [], 1) reshape(p2, [], 1)]);
        corr = R(1,2);
        
        map1(x,y) = corr;
        map2(tx,ty) = max(map2(tx,ty), corr);
    end
end
end
