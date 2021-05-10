clear;
close all;
clc;
rng(4);
kmax = 1;
estimation_accuracy=zeros(kmax,1);
scales = [0.5 0.8 0.9 0.25 0.75 0.85 0.95 1.1 1.2 1.5 1.05 1.15 1.25 1.75 2];
z = size(scales);
for k=1:kmax
    k
    u = randsample(z,1);
    img="Dataset1/im"+k+"/out_s"+scales(u)+"_im"+k+".bmp"
    image = imread(img);
    [r,c,~] = size(image);
%     image = imresize(image,[256,256]);
    %image=rgb2gray(img);
    %image=im2double(img);
    mask = detect(image);
%     mask = imresize(mask,[256,256]);
    mask_im="Dataset1/im"+k+"/out_s"+scales(u)+"_im"+k+"_mask.bmp";
    ref_mask=imread(mask_im);
    imwrite(ref_mask, "generated2/mask_"+k+".bmp");
%     ref_mask = imresize(ref_mask,[256,256]);
    estimation_accuracy(k)= sum((ref_mask-mask).^2, 'all') / (r*c);
end
figure;
semilogy(estimation_accuracy, 'b.-');
figure;
histogram(estimation_accuracy);
