function map = correlation_map2(img,threshold, transform, sz)
[row,col]=size(img);
map=zeros(row,col);
img=padarray(img,[sz,sz],'replicate');

for x=1+sz:row+sz
    for y=1+sz:col+sz
        if map(x-sz,y-sz) ~= 0
            continue;
        end
        
        orig = img(x-sz:x+sz, y-sz:y+sz);
        tran = zeros(1+2*sz,1+2*sz);
        
        oxs = x-sz:x+sz;
        oys = y-sz:y+sz;
        txs = [];
        tys = [];
        
        for dx=-sz:sz
            for dy=-sz:sz
                x_ = round(transform([x+dx;y+dy]));
                tx = x_(1);
                ty = x_(2);
                tx = max(tx,1);
                tx = min(tx,row+2*sz);
                ty = max(ty,1);
                ty = min(ty,col+2*sz);
                tran(dx+sz+1, dy+sz+1) = img(tx,ty);
                txs = [txs tx];
                tys = [tys ty];
            end
        end
        R = corrcoef([reshape(orig, [], 1) reshape(tran, [], 1)]);
        corr = R(1,2);
        
        if(corr>=threshold)
            map(x-sz,y-sz)=1;
            x_ = round(transform([x;y]));
            tx = x_(1);
            ty = x_(2);
            if tx <= sz || tx > row+sz || ty <= sz || ty > col+sz
                continue;
            end
            map(tx-sz,ty-sz)=2;
        end
    end
end
end
