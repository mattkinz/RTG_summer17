function [D,Dt] = FD2D_multiscale(k,levels,p,q)



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
VX = zeros(p,q,levels); VY = VX; 

for ii = 1:levels
    vx = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:q-1)*l(ii)/q)-1).^(k+1))./(exp(-1i*2*pi*(1:q-1)/q)-1)];
    vy = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:p-1)*l(ii)/p)-1).^(k+1))./(exp(-1i*2*pi*(1:p-1)/p)-1)];
    [VX(:,:,ii),VY(:,:,ii)] = meshgrid(vx,vy);
end




D = @(U)D_Forward(U,VX,VY,p,q,levels,k);
Dt = @(Uc)D_Adjoint(Uc,VX,VY,p,q,levels,k);


% multiscale high order finite differences
function [dU] = D_Forward(U,VX,VY,p,q,levels,k)
    
    U = reshape(U,p,q);
    
    % Transform data into frequency domain along each dimension
    % and allocate FD matrices for storage
    Fx = fft(U,q,2); X = zeros(p,q,levels);
    Fy = fft(U,p,1); Y = X;
    
    % filtering for each level and dimension
    for i = 1:levels
        X(:,:,i) = Fx.*VX(:,:,i);
        Y(:,:,i) = Fy.*VY(:,:,i);
    end
    
    % transform back to real space
    X = ifft(X,q,2);
    Y = ifft(Y,p,1);
    
    % reshape data into 3 vectors
    dU = 2^(1-k)/levels*[X(:),Y(:)];
    




% transpose FD
function dtxy = D_Adjoint(dU,VX,VY,p,q,levels,k)
    
    X = reshape(dU(:,1),p,q,levels);
    Y = reshape(dU(:,2),p,q,levels);
    
    % transform data into frequency domain along each dimension
    X = fft(X,q,2);
    Y = fft(Y,p,1);
    
    % conjugate filtering for each level and dimension
    for i = 1:levels
        X(:,:,i) = X(:,:,i).*conj(VX(:,:,i));
        Y(:,:,i) = Y(:,:,i).*conj(VY(:,:,i));
    end
    
    % transform filtered data back to real space
    X = ifft(X,q,2);
    Y = ifft(Y,p,1);
    
    % finish transpose operation by appropriate summing
    X = sum(X,3); Y = sum(Y,3);
    
    dtxy = (X(:) + Y(:))*2^(1-k)/levels;

        

    
    
    
    