function mask=location(map,c)
    [row, col] = size(map);
    
    % filtering image with a 7x7 gaussian o remove noise
    kernel = fspecial('gaussian', [7 7], 1.6);
    mask = imfilter(map, kernel);
    figure;
    imshow(mask);
    idx = (mask>c); % the pixels with value > c are marked 1, others marked 0
    mask(idx) = 1;
    idx = mask~=1;
    mask(idx) = 0;
    
    % removing point clusters having area < 0.1 percent of *(total area)
    area_thresh = floor(0.001*(row*col));
    mask = bwareaopen(mask,area_thresh);
    
    % smoothing out the binary regions  in the mask
    [x,y] = find(mask ==1);
    I1 = zeros(row,col);
    [X,Y] = meshgrid(1:col,1:row) ;
    idx = boundary(x,y) ;
    idx = inpolygon(X(:),Y(:),y(idx),x(idx)) ;
    I1(idx) = 1;
    mask = I1;
    
end
