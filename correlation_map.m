function [mapf,mapb] = correlation_map(img)

[row,col]=size(img);

mapf=zeros(row,col);
mapb=zeros(row,col);

img2=tranform(img);
img3=inv_transform(img);
img1=padarray(img,[2,2],'replicate');
img2=padarray(img2,[2,2],'replicate');
img3=padarray(img3,[2,2],'replicate');

for i=3:row+2
    for j=3:col+2
        patch1=img1(i-2:i+2,j-2:j+2);
        patch2=img2(i-2:i+2,j-2:j+2);
        patch3=img3(i-2:i+2,j-2,j+2);
        
        sum1=sum(patch1.*patch1);
        sum2=sum(patch2.*patch2);
        sum3=sum(patch3.*patch3);
        
        numf=0;
        numb=0;
        
        for i1=1:5
            for j1=1:5
                for i2=1:5
                    for j2=1:5
                        numf=numf+patch1(i1,j1)*patch2(i2,j2);
                        numb=numb+patch1(i1,j1)*patch3(i2,j2);
                    end
                end
            end
        end
        
        mapf(i,j)=numf/(sum1*sum2);
        mapb(i,j)=numb/(sum1*sum3);
        
    end
end

end
