% script for testing accuracy on rotated and copied patches 
clear;
close all;
clc;
rng(1);
kmax = 20;
estimation_accuracy=zeros(kmax,1);
scales = [0.5 0.8 0.9 0.25 0.75 0.85 0.95 1.1 1.2 1.5 1.05 1.15 1.25 1.75 2];
rot = [5 10 15 20 25 30 60 90 120 150 180 210 240 270 300 330];
z = size(scales,2);

for k=1:kmax
    k

    u = randsample(z,1);
    img="Dataset1/im"+k+"/out_r"+rot(u)+"_im"+k+".bmp"
    image = imread(img);
    [r,c,~] = size(image);
    mask = detect(image);
    mask_im="Dataset1/im"+k+"/out_s"+scales(u)+"_im"+k+"_mask.bmp";
    ref_mask=imread(mask_im);
    imwrite(ref_mask, "generated2/mask_"+k+".bmp");
    estimation_accuracy(k)= sum((ref_mask-mask).^2, 'all') / (r*c); % fraction of mismatched points b/w reference mask and the mask produced by this algorithm 
end
figure;
semilogy(estimation_accuracy, 'b.-'); % log scale plots of fraction of mismatched points for each image  
figure; 
histogram(estimation_accuracy);  % histogram of fraction of mismatched points for each image
