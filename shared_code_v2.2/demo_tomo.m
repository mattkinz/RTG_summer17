% demo for simulating tomographic reconstruction using this package
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 05/17/2016


n = 256;  % image dimension
angles = -75:5:75;  % projection angles for tomography

P = phantom(n);  % generate shepp logan phantom
r = radon(P,-angles); %generate radon data using matlab's radon
scale = size(r,1)/n; % scaling needed to generate data matrix
bb = r(:);  %stack the data as a data vector bb
W = radonmatrix(angles,n,size(r,1),scale); % generate sparse tomography matrix


clear pat;  % clear options

% set the desired options into the pat structure for the HOTV code
% for a description of each option, read through check_HOTV_opts
pat.order = 1;
pat.outer_iter = 20;
pat.inner_iter = 10;
pat.mu = 150;
pat.mu0 = 50;
pat.beta = 32;
pat.disp=true;
pat.nonneg=true;
pat.data_mlp = true;
pat.scale_A = true;
pat.scale_b = true;
pat.scale_mu = true;
pat.check_opts = false;
pat.wrap_shrink = false;
pat.L1type = 'isotropic';

% run HOTV code 
[U,~] = HOTV3D(W,bb,[n,n,1],pat);
%% higher order
pat.order = 3;
[U5,~] = HOTV3D(W,bb,[n,n,1],pat);

% multiscale approach
pat.levels = 3;
[U7,~] = HOTV3D(W,bb,[n,n,1],pat);

%% compare with least squares (solved using CGLS)
U2 = run_cgs(W,bb,1e-5,50);
U2 = reshape(U2,n,n);


% compare with algebraic technique using nonnegativity
% (basically just gradient decent)
[U3,~] = SIRT(bb,W,n,100,0);

%% compare with Tikhonov regularized solution
opts.order = 1;
opts.lambda = .25;
opts.scale_lam = true;
opts.maxit = 30;
opts.tol = 1e-3;
ops.scale_A = true;
opts.scale_b = true;
[U4,~] = tikhonov_cgls(W,bb,n,n,1,opts);

% compare with filtered backprojection
U6 = iradon(r,-angles);

%% display results
figure(1);
subplot(2,3,1);imagesc(U,[0 1]);colormap(gray);title('TV');
subplot(2,3,2);imagesc(U5,[0 1]);colormap(gray);title('HOTV3');
subplot(2,3,3);imagesc(U7,[0 1]);colormap(gray);title('MHOTV3');
subplot(2,3,4);imagesc(U3,[0 1]);colormap(gray);title('SIRT');
subplot(2,3,5);imagesc(U4,[0 1]);colormap(gray);title('Tikhonov');
subplot(2,3,6);imagesc(U6,[0 1]);colormap(gray);title('Filerted Backprojection');


figure(2);
subplot(2,1,1);imagesc(P,[0 1]);colormap(gray);title('phantom');
subplot(2,1,2);imagesc(angles,size(r,1),r);colormap(gray);title('sinogram');
xlabel('angle');
