clear;
clc;
estimation_accuracy=zeros(10,1);
for k=1:10
    img="dataset/Dataset_0/im"+k+"_t.bmp";
    image = imread(img);
    image = imresize(image,[256,256]);
    %image=rgb2gray(img);
    %image=im2double(img);
    mask = detect(image);
    mask = imresize(mask,[256,256]);
    mask_im="dataset/Dataset_0/im"+k+"_t_mask.bmp";
    ref_mask=imread(mask_im);
    ref_mask = imresize(ref_mask,[256,256]);
    estimation_accuracy(k)=(sqrt(sum(sum((ref_mask-mask).^2)))/(256*256));
end
figure(1);
plot(estimation_accuracy);
