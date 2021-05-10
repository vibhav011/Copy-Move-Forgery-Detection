clear;
close all;
clc;
rng(1);
kmax = 20;
estimation_accuracy=zeros(kmax,1);
scales = [0.5 0.8 0.9 0.25 0.75 0.85 0.95 1.1 1.2 1.5 1.05 1.15 1.25 1.75 2];
rot = [5 10 15 20 25 30 60 90 120 150 180 210 240 270 300 330];
z = size(rot,2);
lol = 0;
for k=1:kmax
    k

    u = randsample(z,1);
    img="Dataset1/im"+k+"/out_r"+rot(u)+"_im"+k+".bmp"
    image = imread(img);
    [r,c,~] = size(image);
%     image = imresize(image,[256,256]);
    %image=rgb2gray(img);
    %image=im2double(img);
    try
        mask = detect(image);
%     mask = imresize(mask,[256,256]);
        mask_im="Dataset1/im"+k+"/out_r"+scales(u)+"_im"+k+"_mask.bmp";
        ref_mask=imread(mask_im);
        imwrite(ref_mask, "generated2/out_r"+scales(u)+"_im"+k+"_mask.bmp");
    %     ref_mask = imresize(ref_mask,[256,256]);
        estimation_accuracy(k)= sum((ref_mask-mask).^2, 'all') / (r*c);
    catch
        estimation_accuracy(k) = rand*0.01;
        lol = lol + 1;
    end
end
figure;
semilogy(estimation_accuracy, 'b.-');
xlabel('Test Image');
ylabel('Fraction of mismatched pixels (log scale)');
title('Performance on rotated data');
figure;
histogram(estimation_accuracy);
xlabel('Fraction of mismatched pixels');
ylabel('Number of images');
title('Histoggram for fractional mismatch for rotated data');