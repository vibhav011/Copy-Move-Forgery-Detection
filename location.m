function [mask, idx]=location(map, c2)
    [row, col] = size(map);
    
    % lines corresponding to figures can be commented out to monitor progress and detect errors, if any 

    % filtering image with a 7x7 gaussian o remove noise
    kernel = fspecial('gaussian', [7 7], 1.6);
    mask = imfilter(map, kernel);

%     figure;
%     imshow(mask);
    idx = (mask>c2); % the pixels with value > c2 are marked 1, others marked 0
    mask(idx) = 1;
    idx = (mask~=1);
    mask(idx) = 0;
    
%     figure;
%     imshow(mask);
%     drawnow;
    
    % removing point clusters having area < 0.1 percent of *(total area)
    area_thresh = floor(0.001*(row*col));
    mask = bwareaopen(mask,area_thresh);
%     figure;
%     imshow(mask);
    
    % smoothing out the binary regions  in the mask
    [x,y] = find(mask ==1);
    I1 = zeros(row,col);
    [X,Y] = meshgrid(1:col,1:row) ;
    idx2 = boundary(x,y) ; % boundary of masks found 
    idx = inpolygon(X(:),Y(:),y(idx2),x(idx2)) ; % area within boundary filled 
    I1(idx) = 1;
    mask = I1;
%     figure;
%     imshow(mask);
%     drawnow;
    
end
