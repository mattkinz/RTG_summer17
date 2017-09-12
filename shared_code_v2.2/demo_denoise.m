% demo for basic image denoising for RGB image
%
% Written by Toby Sanders @ASU
% School of Math and Stat Sciences
% 06/14/2017

sigma = .1; % noise level

% read image and add noise
X = im2double(imread('skimjamoob13.jpg'));
Xn = X + randn(size(X))*sigma;
[m,n,k] = size(X);

% L1 optimization options
opts.nonneg = true;
opts.mu = 10;
opts.data_mlp = false;
opts.outer_iter = 15;
opts.inner_iter = 10;
opts.tol_inn = 1e-5;
opts.tol_out = 1e-5;
opts.disp = true;
opts.order = 2;

% reconstruct with order 1 (TV) and order 2 (
rec1 = zeros(m,n,k);
rec2 = zeros(m,n,k);
for i = 1:k
    opts.order = 1;
    rec1(:,:,i) = HOTV3D(speye(m*n),col(Xn(:,:,i)),[m,n,1],opts);
    opts.order = 2;
    rec2(:,:,i) = HOTV3D(speye(m*n),col(Xn(:,:,i)),[m,n,1],opts);
end
%% display results
figure(44);
subplot(2,2,1);imagesc(X);title('original')
subplot(2,2,2);imagesc(Xn);title('noisy');
subplot(2,2,3);imagesc(rec1);title('TV denoised');
subplot(2,2,4);imagesc(rec2);title('HOTV denoised');