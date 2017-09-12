%This file is a template demo for running the DART algorithm.
%All of the inputs go into one structure, which is named "opt"
%The options are described at the bottom of this file.
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences 
% 05/2016

n = 256;
angles = -60:5:60;
radius = 30;
n_sigma = 2;


% make the simple phantom image
X = zeros(n);
X(n/2-radius:n/2+radius,n/4:3*n/4)=1;
[i,j] = ind2sub([n,n],1:n^2);
d = (i-(n+1)/2).^2 + (j-(n+1)/2).^2;
s = find(d<radius^2);
X(s)=1/2;


% generate Radon data
r = radon(X,angles);
bb = r(:);
scale = size(r,1)/n;
W = radonmatrix(angles,n,size(r,1),scale); % built radon matrix
noise = random('normal',0,n_sigma,size(bb,1),1);
bb = bb + noise;


clear pat;

pat.order = 1;
pat.outer_iter = 10;
pat.inner_iter = 10;
pat.nonneg = true;
pat.max_c = true;
pat.max_v = 1;
pat.scale_mu = true;
pat.scale_beta = true;
pat.disp = true;


init = HOTV3D(W,bb,[n,n,1],pat);

pat.grays = [0 1/2 1];
pat.init = min(init,1);
pat.mu = 100;
[U,out] = DT_L1(W,bb,[n,n,1],pat);
pat.init = U;
[U2,out2] = DT_L1(W,bb,[n,n,1],pat);

%%
figure(2);
subplot(2,2,1);imagesc(X);title('phantom image');
subplot(2,2,2);imagesc(init);title('TV');
subplot(2,2,3);imagesc(U);title('DT rec');
subplot(2,2,4);imagesc(U2);title('DT rec2');







