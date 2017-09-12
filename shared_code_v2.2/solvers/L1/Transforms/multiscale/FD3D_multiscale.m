function [D,Dt] = FD3D_multiscale(k,levels,p,q,r)


% finite difference operators for higher order TV
% k is the order of the transform
% levels is the number of scales used for the FD transforms
% recommended 3 levels
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/24/2016

l = 0:levels-1;
l = 2.^l;
VX = zeros(p,q,r,levels); VY = VX; VZ = VX;

for ii = 1:levels
    vx = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:q-1)*l(ii)/q)-1).^(k+1))./(exp(-1i*2*pi*(1:q-1)/q)-1)];
    vy = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:p-1)*l(ii)/p)-1).^(k+1))./(exp(-1i*2*pi*(1:p-1)/p)-1)];
    vz = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:r-1)*l(ii)/r)-1).^(k+1))./(exp(-1i*2*pi*(1:r-1)/r)-1)];
    [VX(:,:,:,ii),VY(:,:,:,ii),VZ(:,:,:,ii)] = meshgrid(vx,vy,vz);
end


D = @(U)D_Forward(U,VX,VY,VZ,p,q,r,levels,k);
Dt = @(Uc)D_Adjoint(Uc,VX,VY,VZ,p,q,r,levels,k);


% multiscale high order finite differences
function [dU] = D_Forward(U,VX,VY,VZ,p,q,r,levels,k)

    U = reshape(U,p,q,r);
    
    % Transform data into frequency domain along each dimension
    % and allocate FD matrices for storage
    Fx = fft(U,q,2); X = zeros(p,q,r,levels);
    Fy = fft(U,p,1); Y = X;
    Fz = fft(U,r,3); Z = X;
    
    % filtering for each level and dimension
    for i = 1:levels
        X(:,:,:,i) = Fx.*VX(:,:,:,i);
        Y(:,:,:,i) = Fy.*VY(:,:,:,i);
        Z(:,:,:,i) = Fz.*VZ(:,:,:,i);
    end
    
    % transform back to real space
    X = ifft(X,q,2);
    Y = ifft(Y,p,1);
    Z = ifft(Z,r,3);
    
    % reshape data into 3 vectors
    dU = 2^(1-k)/levels*[X(:),Y(:),Z(:)];
    




% transpose FD
function dtxy = D_Adjoint(dU,VX,VY,VZ,p,q,r,levels,k)
    
    X = reshape(dU(:,1),p,q,r,levels);
    Y = reshape(dU(:,2),p,q,r,levels);
    Z = reshape(dU(:,3),p,q,r,levels);
    
    % transform data into frequency domain along each dimension
    X = fft(X,q,2);
    Y = fft(Y,p,1);
    Z = fft(Z,r,3);
    
    % conjugate filtering for each level and dimension
    for i = 1:levels
        X(:,:,:,i) = X(:,:,:,i).*conj(VX(:,:,:,i));
        Y(:,:,:,i) = Y(:,:,:,i).*conj(VY(:,:,:,i));
        Z(:,:,:,i) = Z(:,:,:,i).*conj(VZ(:,:,:,i));
    end
    
    % transform filtered data back to real space
    X = ifft(X,q,2);
    Y = ifft(Y,p,1);
    Z = ifft(Z,r,3);
    
    % finish transpose operation by appropriate summing
    X = sum(X,4); Y = sum(Y,4); Z = sum(Z,4);
    
    dtxy = 2^(1-k)/levels*(X(:) + Y(:) + Z(:));

        

    
    
    
    