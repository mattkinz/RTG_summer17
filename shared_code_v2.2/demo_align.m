% demo to try out the alignment algorithms for aligning electron tomography
% data.  Main work is the center of mass method found here: 
% http://ascimaging.springeropen.com/articles/10.1186/s40679-015-0005-7
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 06/08/2016



d = 512;  % image dimension
angles = -80:8:80;  % angular sampling
noise_mu = 0;   % mean noise level
noise_sigma = 5;  % standard deviation


% generate the phantom and tomography data, and add in noise
P = phantom(d);   
R = radon(P,-angles);
scale = size(R,1)/d;  

noise = random('normal',noise_mu,noise_sigma,size(R,1),size(R,2));
R = R + noise;

% display noisy sinogram
figure(1);subplot(2,2,1);
imagesc(angles,size(R,1),R);title('Noisy sinogram');
xlabel('angle ');colormap(gray);

% reshape the data so that the angle information changes along the 3rd
% dimension
R = reshape(R,size(R,1),1,size(R,2));


% generate random shifts for simulate misaligned data
rshift = 20*rand(size(R,3),1);
rshift = round(rshift);
for i = 1:size(R,3)
    R(:,1,i) = circshift(R(:,1,i),[rshift(i),0]);
end

% display misaligned data
figure(1);subplot(2,2,2);
imagesc(angles,size(R,1),squeeze(R));title('Noisy, misaligned sinogram');
xlabel('angle');colormap(gray);



% apply center of mass and cross correlation alignment methods
[Rcom,scom] = COM_align(R,angles,1/2,size(R,3));
[Rcc,scc] = cross_corr(R);

figure(1);subplot(2,2,3);
imagesc(angles,size(R,1),squeeze(Rcom));title('COM aligned sinogram');
xlabel('angle');colormap(gray);
subplot(2,2,4);
imagesc(angles,size(R,1),squeeze(Rcc));title('cross corr aligned signogram');
xlabel('angle');colormap(gray);
% compute reconstructions for each alignment and compare results
% SIRT reconstructions
W = radonmatrix(angles,d,size(R,1),scale);
rec_com = SIRT(Rcom,W,d,50,0);
rec_unali = SIRT(R,W,d,50,0);
rec_cc = SIRT(Rcc,W,d,50,0);
figure(2);
subplot(2,2,1);imagesc(rec_com,[0 1]);title('reconstruction, COM alignment');colormap(gray);
subplot(2,2,2);imagesc(rec_unali,[0 1]);title('reconstruction, no alignment');colormap(gray);
subplot(2,2,3);imagesc(rec_cc,[0 1]);title('reconstruction, cross corr alignment');colormap(gray);
subplot(2,2,4);imagesc(P,[0 1]);title('phantom');colormap(gray);