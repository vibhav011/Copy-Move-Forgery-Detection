function [map] = correlation_alternate(img,threshold)

[row,col]=size(img);

map=zeros(row,col);

img=padarray(img,[2,2],'replicate');

for x=3:row+2
    for y=3:col+2 
        [tx,ty]=transform(x,y);
        num=0;
        den1=0;
        den2=0;
        for dx=-2:2
            for dy=-2:2
              den1=den1+img(x+dx,y+dy)*img(x+dx,y+dy);
              den2=den2+img(tx+dx,ty+dy)*img(tx+dx,ty+dy);
                for tdx=-2:2
                    for tdy=-2:2
                      num=num+img(x+dx,y+dy)*img(tx+tdx,ty+tdy);
                    end
                end
            end
        end
        den=den1*den2;
        corr=num/den;
        if(corr>=threshold)
            map(x,y)=1;
            map(tx,ty)=1;
        end
    end
end
end

end
