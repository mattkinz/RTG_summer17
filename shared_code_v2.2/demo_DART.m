%This file is a template demo for running the DART algorithm.
%All of the inputs go into one structure, which is named "opt"
%The options are described at the bottom of this file.
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 05/2016

n = 256;
angles = -70:2:70;
radius = 30;


% make the phantom image
X = zeros(n);
X(n/2-radius:n/2+radius,n/4:3*n/4)=1;
[i,j] = ind2sub([n,n],1:n^2);
d = (i-(n+1)/2).^2 + (j-(n+1)/2).^2;
s = find(d<radius^2);
X(s)=1/2;
%imagesc(X);pause;

% generate Radon data
r = radon(X,angles);
bb = r(:);
scale = size(r,1)/n;
W = radonmatrix(angles,n,size(r,1),scale); % built radon matrix
Uinit = SIRT(bb,W,n,50,0);  % compute solution with SIRT

clear opt;


%%%% MAIN INPUTS %%%%%%
opt.angles=angles;
opt.thresh=[.25 .75];
opt.grays=[0 .5 1];
%opt.data_name = 'stack-ali.mat';
%opt.data_type='mat';

opt.resolution=n;
opt.startslice=1;
opt.endslice=1;



% set iterations for each refinement
opt.bdry_dirt_iter=30;
opt.region_dirt_iter=0;
opt.t_tol=.01;
opt.t_delta=.005;
opt.t_epsilon=.01;
opt.bdry_update_type='SIRT';
opt.region_update_type='SIRT';
opt.rand_pixel_probability = 0;

%If the an initial solution exists
opt.initialsolution=false;
opt.initialfile='sirt-rec.mat';



%%%  ADDITIONAL USER OPTIONS
opt.W=W;
opt.model='line';
opt.inner_iter=10;
opt.minc=true;
opt.disp=true;
opt.disppic=true;
opt.convergence_criteria=false;
opt.convergence_tol=10^(-6);
opt.max_iter=60;
opt.rx=5;
opt.rz=5;
opt.sigmax=3;
opt.sigmaz=1.5;
opt.chunksize=50;
opt.overlap=10;





%run DART
[Ud,out,~]=DIRT(opt,bb,Uinit);

%figure(2);imagesc(X);title('phantom image');
%figure(3);imagesc(Ud);title('DART reconstruction');
%% display results
figure;
subplot(2,2,1);imagesc(X);title('phantom image');colormap(gray);
subplot(2,2,2);imagesc(Uinit,[0 1]);title('SIRT');colormap(gray);
subplot(2,2,3);imagesc(threshold(Uinit,[.25 .75],[0 1/2 1]));colormap(gray);
title('SIRT segmented');
subplot(2,2,4);imagesc(Ud);title('DART');colormap(gray);
%SAVE THE RECONSTRUCTION
%save('dart-stuff','dartrec','out','initialrec','-v7.3');









%OUTPUTS:
    %U - the 3-D DART reconstruction.  The slices change with the 3rd
        %dimension, i.e. dartrec(:,:,i) holds the ith slice of the
        %reconstruction
    %out - data about the reconstruction
    %init - the SIRT reconstruction or input reconstruction used as
        %the initial solution to compute "dartrec".


%INPUTS:
%Inputs that must be user defined:
    %opt.angles - the projection angles listed in order, in degrees
    %opt.thresh - a vector holding the thresholds used for segmenting 
    %opt.grays - a vector holding the gray values used in segmentation.
    %opt.data_name - name of the file in which the aligned tilt series is 
        %saved.  Note the tilt axis should be horizontal and centered.
    %opt.data_type - type of file of opt.filename, set to either 'mrc' or 
        %'mat'

%Inputs recommended to be user defined, but default values are otherwise
%used:
    %opt.bdry_dirt_iter - number of DART iterations with boundary cleaning
    %opt.region_dirt_iter - number of region refinements (DIP-LS)
    %opt.t_tol - radius of interval used for the partial segmentation
        %function with DIPS.  Set a small tolerance for slower more
        %careful convergence and a high tolerance for fast convergence.
        %Default is set accordingly to the gray levels.
    %opt.t_delta - increase in opt.t_tol whenever new pixels are not being
        %classified.  Default is opt.t_tol*1/4;
    %opt.t_epsilon - this epsilon is used to determine if opt.t_tol should 
        %be increased by opt.t_delta.  It is the case if the number of 
        %classified pixels has only changed relatively by less that
        %epsilon.  Default is .01.
    %opt.bdry_update_type - indicates the solver used for bdry updates.
        %Set to 'SIRT' (default) or 'CGLS'.
    %opt.region_update_type - indicates the solver used for region updates.
        %Set to 'SIRT' (default) or 'CGLS'. 
    %opt.rand_pixel_probability - probability to a free pixel is set to
        %free.  Default is 0.
    %opt.resolution - the number of pixels used for reconstruction.
        %Reconstruction slices will be opt.recsize by opt.recsize. Default
        %value is the detector count
    %opt.startslice - which slice the reconstruction begins on.  Default is 1. 
    %opt.endslice - which slice the reconstruction ends on.  Default is the
    %end of the stack.
    %opt.initialsolution - Set to "true" an initial solution has been
        %computed and will be used, otherwise set to "false".
    %opt.initialfile - name of the matlab file where the reconstruction is
        %saved
        
%More advanced inputs, and default values are again used these are not
%specified:
    %opt.inner_iter - number of SIRT iterations in each DART iteration.
        %Default value is 10.
    %opt.minc - option to use the density constraint within the solver.
        %Set to true or false. Default value is true.
    %opt.W - if the projection matrix is precomputed, set opt.W to be
        %this projection matrix.  Otherwise set to false.
    %opt.disp - prints information about the reconstrution at each
        %iteration.  Set to "true" for the display and "false" 
        %otherwise.  Default is "true".
    %opt.disppic - displays various images about the reconstruction at
        %each iteration. Set to "true" for this display and "false"
        %otherwise.  Default is "true".
    %opt.convergence_criteria - if set to "true", DART will iterate until
        %the reconstruction converges.  Set to "false" to simply use
        %opt.outeriter instead.  Default is "false".  Typically
        %convergence is seen after just a few iterations.
    %opt.convergence_tol - tolerance for opt.convergencecriteria.
        %Default is 10^(-3).
    %opt.max_iter - maximum number of iterations if
        %opt.convergencecriteria is used.  Default is 50.
    %opt.rx - radius of the gaussian in the slice (x,y) dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.rz - radius of the gaussian in the slice z dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.sigmax - sigma value used for the gaussian smoothing in the
        %(x,y) dimension.  Default is 1.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.sigmaz - sigma value used for the gaussian smoothing in the
        %z dimension.  Default is 1.5.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.chunksize - size of the chunks that DART will run on.  Default
        %is 50.
    %opt.overlap - once the chunks are computed, they are merged
        %together using a partition of unity.  The overlap is the overlap
        %of this merging.  Default is 10.



        
        
%%%%%%Suggested default options%%%%%%
%{
opt.inneriter=10;
opt.outeriter=10;
opt.convergencecriteria=1;
opt.convergencetol=10^(-3);
opt.maxiter=60;
opt.rx=5;
opt.rz=5;
opt.sigmax=2;
opt.sigmaz=1.5;
opt.chunksize=50;
opt.overlap=10;
opt.regionrefinements=1;
opt.includebdry=1;
opt.inittol=.2;
opt.finaltol=.2;
opt.whenstopregion=25;
%}


