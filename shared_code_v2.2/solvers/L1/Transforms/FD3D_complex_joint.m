function [D,Dt] = FD3D_complex_joint(kx,kz,p,q,r,theta)



% finite difference operators for higher order TV
% k is the order of the transform
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

D = @(U)D_Forward(U,kx,kz,p,q,r,theta);
Dt = @(Uc)D_Adjoint(Uc,kx,kz,p,q,r,theta);


% high order finite differences
function dU = D_Forward(U,kx,kz,p,q,r,theta)
U = reshape(U,p,q,r);
U = U.*theta;
dU = zeros(p,q,r,3);

if kx~=0
    if kx<=q
        dU(:,:,:,1) = diff([U,U(:,1:kx,:)],kx,2);
    end
    if kx<=p
        dU(:,:,:,2) = diff([U;U(1:kx,:,:)],kx,1);
    end
else
    %standard l1 minimization of order 0
    dU(:,:,:,1) = U;
    dU(:,:,:,2) = U;
end
if kz~=0
    if kz<=r
        dU(:,:,:,3) = diff(cat(3,U,U(:,:,1:kz)),kz,3);
    end
else
    dU(:,:,:,3) = U;
end
dU = reshape(dU,p*q*r,3);
    



%transpose FD
function U = D_Adjoint(dU,kx,kz,p,q,r,theta)
U = zeros(p,q,r);
dU = reshape(dU,p,q,r,3);
    
if kx~=0 
    if kx<=q
        U = U + (-1)^kx*diff([dU(:,end-kx+1:end,:,1),dU(:,:,:,1)],kx,2);
    end    
    if kx<=p
        U = U + (-1)^kx*diff([dU(end-kx+1:end,:,:,2);dU(:,:,:,2)],kx,1);
    end 
else
    %standard l1 minimization
    U = sum(dU(:,:,:,1:2),4);
end

if kz~=0
    if kz<=r
        U = U + (-1)^kz*diff(cat(3,dU(:,:,end-kz+1:end,3),dU(:,:,:,3)),kz,3);
    end
else
    U = U + dU(:,:,:,3);
end

U = U.*conj(theta);
U = U(:);
    
    