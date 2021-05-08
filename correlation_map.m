function map = correlation_map(img,threshold, transform)
[row,col]=size(img);
map=zeros(row,col);
img=padarray(img,[2,2],'replicate');

for x=3:row+2
    for y=3:col+2
        if map(x-2,y-2) ~= 0
            continue;
        end
        num=0;
        den1=0;
        den2=0;
        for dx=-2:2
            for dy=-2:2
                x_ = round(transform([x+dx;y+dy])); % threshold returns the affine transoformed pixel location as a 2-tuple
                tx = x_(1);
                ty = x_(2);
                if tx <= 0 || tx > row+4 || ty <= 0 || ty > col+4
                    continue;
                end
                num=num+img(x+dx,y+dy)*img(tx,ty);
                den1=den1+img(x+dx,y+dy)^2;
                den2=den2+img(tx,ty)^2;
            end
        end
        den=sqrt(double(den1*den2));
        corr=num/den;
        if(corr>=threshold)
            map(x-2,y-2)=1;
            x_ = round(transform([x;y]));
            tx = x_(1);
            ty = x_(2);
            if tx <= 2 || tx > row+2 || ty <= 2 || ty > col+2
                continue;
            end
            map(tx-2,ty-2)=2;
        end
    end
end
end
