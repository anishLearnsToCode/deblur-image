% function [H] = kurt_bld_dconv(Ib,kernel_size,sample_rate,alpha,lambda)
clc;
clear;
close all;


image_name = "bridge";
image_format = "png";
source_destination = "../input/" + image_name + "." + image_format;
I = imread(source_destination);

[YR, kernel_r] = deblur_single_channel(image_name, image_format, 1);
[YG, kernel_g] = deblur_single_channel(image_name, image_format, 2);
[YB, kernel_b] = deblur_single_channel(image_name, image_format, 3);
N = cat(3, YR, YG, YB);

figure(10);
subplot(121); imshow(I);
subplot(122); imshow(N);
