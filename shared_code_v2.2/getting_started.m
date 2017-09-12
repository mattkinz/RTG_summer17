% This package includes algorithms primarily for solving inverse problems 
% with a particular focus on tomographic reconstruction and regularization.
% Also included are image alignment algorithms
%
% To get started, please look at the read-me PDF file, AND run this script 
% to add these folders to your MATLAB paths so that everything works
% Demos are included for each code of interest
%
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 12/15/2016
%
% Please email any bugs or inquires to toby.sanders@asu.edu

mainpath = pwd;

addpath([mainpath,'/solvers']);
addpath([mainpath,'/solvers/L1']);
addpath([mainpath,'/solvers/L1/Transforms']);
addpath([mainpath,'/solvers/L1/Transforms/multiscale']);
addpath([mainpath,'/solvers/L1/Transforms/wave_shear']);
addpath([mainpath,'/solvers/L1/utilities']);
addpath([mainpath,'/solvers/L1/inpaint']);
addpath([mainpath,'/solvers/tikhonov']);
addpath([mainpath,'/tomography/SIRT']);
addpath([mainpath,'/tomography/JOHANN_RADON']);
addpath([mainpath,'/tomography/Align']);
addpath([mainpath,'/tomography/Align/subroutines'])
addpath([mainpath,'/tomography/DART']);
addpath([mainpath,'/tomography/DART/subroutines'])
addpath([mainpath,'/utilities']);