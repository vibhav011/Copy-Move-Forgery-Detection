function [map] = correlation_map(img,threshold)
[row,col]=size(img);
map=zeros(row,col);
img=padarray(img,[2,2],'replicate');
for x=3:row+2
    for y=3:col+2  
        num=0;
        den1=0;
        den2=0;
        for dx=-2:2
            for dy=-2:2
                [tx,ty]=transform(x+dx,y+dy); % threshold returns the affine transoformed pixel location as a 2-tuple
                num=num+img(x+dx,y+dy)*img(tx,ty);
                den1=den1+img(x+dx,y+dy)*img(x+dx,y+dy);
                den2=den2+img(tx,ty)*img(tx,ty);
            end
        end
        den=sqrt(den1*den2);
        corr=num/den;
        if(corr>=threshold)
            map(x,y)=1;
            [tx,ty]=transform(x,y);
            map(tx,ty)=1;
        end
    end
end
end
