% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 11/2016

d = 400;  % image dimension
subsamp =.2; % subsampling rate

P = phantom(d);  % generate phantom


% set l1 optimization parameters 
pat.mu = 200;
pat.inner_iter = 15;
pat.outer_iter = 30;
pat.nonneg=true;
pat.disp = true;
pat.mu0 = 50;
pat.order = 1;
pat.levels = 1;
pat.data_mlp = true;
pat.tol_inn = 1e-5;
pat.tol_out = 1e-5;
pat.L1type = 'isotropic';


% generate random samples
scnt = round(subsamp*d^2);
S = rand(d^2,1);
[~,S] = sort(S);
S = S(1:scnt);  % sample indices
bb = P(S);  % known samples
X = zeros(size(P));
X(S) = P(S);

%% run HOTV inpainting algorithm
[U,out] = inpaint_3D(bb,S,[d,d,1],pat);
figure;
subplot(2,2,1);imagesc(P);title('original');
subplot(2,2,2);imagesc(X);title('subsampled');
subplot(2,2,3);imagesc(U);title('inpainted');