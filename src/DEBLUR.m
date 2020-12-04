% function [H] = kurt_bld_dconv(Ib,kernel_size,sample_rate,alpha,lambda)
clc;
clear;
close all;

kernel_size = [3 3];
sample_rate = 0.75; 
lambda = 0.002;
alpha = 5500;

image_name = "lenna";
image_format = "png";
source_destination = "../input/" + image_name + "." + image_format;
target_destination = get_target_destination(image_name, image_format);
T = imread(target_destination);
I = imread(source_destination);
y = I;
B = im2double(rgb2gray(I));
Ib = B;


figure(1)
subplot(221)
imshow(I)

%%
[rsize_I,csize_I] = size(Ib);
rsize_ker = kernel_size(1);
csize_ker = kernel_size(2);
len_ker = rsize_ker*csize_ker;
rsize_sample = 2*round(sample_rate*rsize_ker) + 1;
csize_sample = 2*round(sample_rate*csize_ker) + 1;
rrad_ker = round((rsize_ker-1)/2);
crad_ker = round((csize_ker-1)/2);
saveres = false;
% if saveres
%     sdir = ['temp/' num2str(round(rand()*100000000))];
%     mkdir(sdir);
% end
%% 
% disp('computing the Hessian matrix ...');
%extract features 
g = fspecial('unsharp');
B = conv2(Ib,g,'same');
g = fspecial('log');
B = conv2(B,g,'same');


% 'kurtosis'
[K,patch,v]=kurtosis_var_patch(B);

H=0;
for i=1:size(patch,1)
    
    A_p=(convmtx(patch(i,:),kernel_size(1)*kernel_size(1)));
    
    Ai=power(A_p,1);
    Ai=Ai/max(Ai(:));
%      Ai = reshape(Ai,rsize_sample,csize_sample);
    

%     Ai = conv2mtx(Ai,rsize_ker,csize_ker,'same');
    H = H + power((Ai*Ai')/(K(i)*K(i)*v(i)*1000000),1);

end
    H=power(H,1);
    
I0 = Ib;

rm1 = min(2*rsize_ker,round(rsize_I/5));
cm1 = min(2*csize_ker,round(csize_I/5));
rm = rm1+rrad_ker;
cm = cm1+crad_ker;
vB = Ib(rm + 1:rsize_I - rm, cm + 1:csize_I - cm);
vB = reshape(vB,size(vB,1)*size(vB,2),1);

options = optimset('quadprog');
options.Diagnostics = 'off';
options.Display = 'off';
options.LargeScale = 'off';
options.Algorithm = 'interior-point-convex';
Aeq = ones(1,len_ker);
beq = 1;

lb = zeros(len_ker,1);
ub = ones(len_ker,1);
x0 = ub./sum(ub);
disp('loop ...');
iter = 0;
lambda = 1/lambda;

for level=1:1
    if level == 1
        max_iter = 10;
        lambda1 = lambda/20;
    elseif level == 2
        max_iter = 20;
        lambda1 = lambda/15;
    elseif level ==3 
        max_iter = 35;
        lambda1 = lambda/10;
    elseif level == 4
        max_iter = 50;
        lambda1 = lambda/7;
    elseif level ==5 
        max_iter = 70;
        lambda1 = lambda/3;
    else
        max_iter = 200;
        lambda1 = lambda;
    end
    obj0 = inf;
    while iter < max_iter
        iter = iter + 1;

        I01 = I0(rm1+1:rsize_I-rm1, cm1+1:csize_I-cm1);
        
         A = conv2mtx(I01,rsize_ker,csize_ker,'valid');
         A=power(A,1);
         
        Ai=Ai/max(Ai(:));
         ku=kurtosis(I01(:));
        va=var(I01(:));
        f =  - A'*vB;
%         G = alpha*H + power(A'*A/(ku*ku*va*va*100000),1);
G = alpha*H + power(A'*A,1);
        [z,obj] = quadprog(G,f,[],[],Aeq,beq,lb,ub,x0,options);

        kernel = reshape(z,rsize_ker,csize_ker);
%         kernel=kernel/max(kernel(:));
%         kernel=1-kernel;
        kernel=kernel/sum(kernel(:));
 
        kernel=centerPSF(kernel,40/255);
        
                I0 = fast_deconv(Ib, kernel, lambda1, 1, I0);
%                 I0 = deconvblind(Ib,kernel);
  obj = obj + 0.5*norm(vB).^2;

figure(1)
subplot(222)
imshow(kernel,[])

%         disp([' iter ' num2str(iter) ', obj=' num2str(obj)]);
         disp([' iter ' num2str(iter)]);
        if abs(obj - obj0) < 1e-6*obj 
            break;
        else
            obj0 = obj;
        end        
        if saveres && mod(iter,5) == 1
            imwrite(kernel./max(kernel(:)),[sdir '/iter' num2str(iter) '_ker.png']);
        end
    end
end

kernel=kernel/max(kernel(:));
% kernel=1-kernel;
% 
% x=kernel;
% x=x/max(x(:));
% for i=1:13
%     for j=1:13
%         if x(i,j)>.7
%             x(i,j)=(x(i,j)*x(i,j));
%         else
%             x(i,j)=x(i,j)/10;
%         end
%     end
% end
% kernel=x;
%         kernel=kernel/max(kernel(:));
% 
kernel=kernel/sum(kernel(:));
 parameters;
yout = fftCGSRaL(Ib,kernel,PAR);
% % yout = fast_deconv(B, k, lambda, 1, B);
% 
% yout=uint8(yout);
figure(1)
subplot(221)
imshow(I)
subplot(222)
imshow(T,[])
subplot(223)
imshow(kernel,[])
subplot(224)
imshow(yout,[])