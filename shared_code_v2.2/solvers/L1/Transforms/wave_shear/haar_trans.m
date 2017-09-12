function [D,Dt] = haar_trans(p,q)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


Z = zeros(p,q);Z(1,1) = 1;Z(2,2)=1;Z(1,2) = -1;Z(2,1)=-1;
Z = fft2(Z);

D = @(U)haar_forward(U,Z,p,q);
Dt =@(dU)haar_adj(dU,Z,p,q);

function U = haar_forward(U,Z,p,q)

U = fft2(reshape(U,p,q));
U = (ifft2(U.*conj(Z)));


function dU = haar_adj(dU,Z,p,q)

dU = col(ifft2(fft2(dU).*Z));


