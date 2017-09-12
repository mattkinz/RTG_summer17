function [D,Dt] = FD3D_weighted(k,p,q,r,c)



% finite difference operators for higher order TV
% k is the order of the transform
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

D = @(U)D_Forward(U,k,p,q,r,c);
Dt = @(Uc)D_Adjoint(Uc,k,p,q,r,c);


% high order finite differences
function [dU] = D_Forward(U,k,p,q,r,c)
if k~=0
    U = reshape(U,p,q,r);
    dU = zeros(p,q,r,3);
    if k<=q
        dU(:,:,:,1) = diff([U,U(:,1:k,:)],k,2);
    end
    if k<=p
        dU(:,:,:,2) = diff([U;U(1:k,:,:)],k,1);
    end
    if k<=r
        dU(:,:,:,3) = diff(cat(3,U,U(:,:,1:k)),k,3);
    end
    dU = c.*reshape(dU,p*q*r,3);
else
    %standard l1 minimization of order 0
    dU = c.*U(:);
end

dU = dU*2^(1-k);  % normalization

%transpose FD
function U = D_Adjoint(dU,k,p,q,r,c)

dU = dU.*conj(c);

if k~=0
    U = zeros(p,q,r);
    dU = reshape(dU,p,q,r,3);
    if k<=q
        U = U + (-1)^k*diff([dU(:,end-k+1:end,:,1),dU(:,:,:,1)],k,2);
    end    
    if k<=p
        U = U + (-1)^k*diff([dU(end-k+1:end,:,:,2);dU(:,:,:,2)],k,1);
    end
    if k<=r
        U = U + (-1)^k*diff(cat(3,dU(:,:,end-k+1:end,3),dU(:,:,:,3)),k,3);
    end
    U = U(:);
else
    %standard l1 minimization
    U = dU(:);
end
    
    
U = U*2^(1-k);  % normalization
    
    
    