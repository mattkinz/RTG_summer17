function [D,Dt] = FFD2D_multiscale(k,levels,p,q)


% finite difference operators for higher order TV
% k is the order of the transform
% levels is the number of scales used for the FD transforms
% recommended 3 levels
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/24/2016


D = @(U)D_Forward(U,k,p,q,levels);
Dt = @(dU)D_Adjoint(dU,k,p,q,levels);

function dU = D_Forward(U,k,p,q,levels)


kq = min(k,q);
kp = min(k,p);

U = reshape(U,p,q);
dU = zeros(p,q,levels,2);
%D^k*U, in 2-D
dU(:,:,1,1) = diff([U,U(:,1:kq)],kq,2);
dU(:,:,1,2) = diff([U;U(1:kp,:)],kp,1);


% upscaling levels for lower frequencies
% P^(k+1)*(D^k)*U, in 3-D
for i = 2:levels
    for j = 0:k+1
        s = j*(i-1);
        if j<=kq+1
            dU(:,:,i,1) = dU(:,:,i,1) ...
                + nchoosek(kq+1,j)*circshift(dU(:,:,i-1,1),[0 -s 0]);
        end
        if j<=kp+1
            dU(:,:,i,2) = dU(:,:,i,2) ...
                + nchoosek(kp+1,j)*circshift(dU(:,:,i-1,2),[-s 0 0]);
        end       
    end
    dU(:,:,i,:) = 0.5*dU(:,:,i,:);
end
dU = dU(:)*2^(1-k)/levels;

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjoint operation

function U = D_Adjoint(dU,k,p,q,levels)

kq = min(k,q);
kp = min(k,p);
dU = reshape(dU,p,q,levels,2);
U = reshape(dU(:,:,1,:),p,q,2);
dU(:,:,1,:)='';
for i = 2:levels
    dUp = zeros(size(dU));
    for j = 0:k+1
        s = j*(i-1);
        if j<=kq+1
           dUp(:,:,:,1) = dUp(:,:,:,1) ...
               + nchoosek(kq+1,j)*circshift(dU(:,:,:,1),[0 s 0]);
        end
        if j<=kp+1
           dUp(:,:,:,2) = dUp(:,:,:,2) ...
               + nchoosek(kp+1,j)*circshift(dU(:,:,:,2),[s 0 0]);
        end     
    end
    U(:,:,1) = U(:,:,1) + 2^(1-i)*dUp(:,:,1,1);
    U(:,:,2) = U(:,:,2) + 2^(1-i)*dUp(:,:,1,2);
    dU = dUp(:,:,2:end,:);
end

U(:,:,1) = (-1)^kq*diff([U(:,end-kq+1:end,1),U(:,:,1)],kq,2);
U(:,:,2) = (-1)^kp*diff([U(end-kp+1:end,:,2);U(:,:,2)],kp,1);


U = col(sum(U,3))*2^(1-k)/levels;



