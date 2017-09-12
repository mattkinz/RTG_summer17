function [D,Dt] = FFD3D_multiscale(k,levels,p,q,r)


% finite difference operators for higher order TV
% k is the order of the transform
% levels is the number of scales used for the FD transforms
% recommended 3 levels
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/24/2016


D = @(U)D_Forward(U,k,p,q,r,levels);
Dt = @(dU)D_Adjoint(dU,k,p,q,r,levels);

function dU = D_Forward(U,k,p,q,r,levels)


kq = min(k,q);
kp = min(k,p);
kr = min(k,r);

U = reshape(U,p,q,r);
dU = zeros(p,q,r,levels,3);
%D^k*U, in 3-D
dU(:,:,:,1,1) = diff([U,U(:,1:kq,:)],kq,2);
dU(:,:,:,1,2) = diff([U;U(1:kp,:,:)],kp,1);
dU(:,:,:,1,3) = diff(cat(3,U,U(:,:,1:kr)),kr,3);


% upscaling levels for lower frequencies
% P^(k+1)*(D^k)*U, in 3-D
for i = 2:levels
    for j = 0:k+1
        s = j*(i-1);
        if j<=kq+1
            dU(:,:,:,i,1) = dU(:,:,:,i,1) ...
                + nchoosek(kq+1,j)*circshift(dU(:,:,:,i-1,1),[0 -s 0]);
        end
        if j<=kp+1
            dU(:,:,:,i,2) = dU(:,:,:,i,2) ...
                + nchoosek(kp+1,j)*circshift(dU(:,:,:,i-1,2),[-s 0 0]);
        end       
        if j<=kr+1
            dU(:,:,:,i,3) = dU(:,:,:,i,3) ...
                + nchoosek(kr+1,j)*circshift(dU(:,:,:,i-1,3),[0 0 -s]);
        end
    end
    dU(:,:,:,i,:) = 0.5*dU(:,:,:,i,:);
end
dU = 2^(1-k)/levels*dU(:);
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjoint operation

function U = D_Adjoint(dU,k,p,q,r,levels)

kq = min(k,q);
kp = min(k,p);
kr = min(k,r);
dU = reshape(dU,p,q,r,levels,3);
U = reshape(dU(:,:,:,1,:),p,q,r,3);
dU(:,:,:,1,:)='';
for i = 2:levels
    dUp = zeros(size(dU));
    for j = 0:k+1
        s = j*(i-1);
        if j<=kq+1
           dUp(:,:,:,:,1) = dUp(:,:,:,:,1) ...
               + nchoosek(kq+1,j)*circshift(dU(:,:,:,:,1),[0 s 0]);
        end
        if j<=kp+1
           dUp(:,:,:,:,2) = dUp(:,:,:,:,2) ...
               + nchoosek(kp+1,j)*circshift(dU(:,:,:,:,2),[s 0 0]);
        end
        if j<=kr+1
           dUp(:,:,:,:,3) = dUp(:,:,:,:,3) ...
               + nchoosek(kr+1,j)*circshift(dU(:,:,:,:,3),[0 0 s]);
        end        
    end
    U(:,:,:,1) = U(:,:,:,1) + 2^(1-i)*dUp(:,:,:,1,1);
    U(:,:,:,2) = U(:,:,:,2) + 2^(1-i)*dUp(:,:,:,1,2);
    U(:,:,:,3) = U(:,:,:,3) + 2^(1-i)*dUp(:,:,:,1,3);
    dU = dUp(:,:,:,2:end,:);
end

U(:,:,:,1) = (-1)^kq*diff([U(:,end-kq+1:end,:,1),U(:,:,:,1)],kq,2);
U(:,:,:,2) = (-1)^kp*diff([U(end-kp+1:end,:,:,2);U(:,:,:,2)],kp,1);
U(:,:,:,3) = (-1)^kr*diff(cat(3,U(:,:,end-kr+1:end,3),U(:,:,:,3)),kr,3);


U = 2^(1-k)/levels*col(sum(U,4));


