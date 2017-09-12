function [U,out] = denoise_3D(I,opts)


% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/22/2016


A = @(x,mode)Idenity_map(x,mode);
[p,q,r] = size(I);
[U,out] = PA3D(A,I(:),[p,q,r],opts);



function x = Idenity_map(x,mode)

x = x(:);