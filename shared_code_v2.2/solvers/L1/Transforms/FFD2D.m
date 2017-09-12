function [D,Dt] = FFD2D(k,p,q,r)


% finite difference operators for higher order TV
% k is the order of the transform
%
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016



VX = (exp(1i*2*pi*(0:q-1)/q)-1).^k;
VY = (exp(1i*2*pi*(0:p-1)/p)-1).^k;
VZ = (exp(1i*2*pi*(0:r-1)/r)-1).^k;
%VX = ([1:q/2,q/2:-1:1]/(q/2)).^k;
%VY = ([1:p/2,p/2:-1:1]/(p/2)).^k;
%VZ = [1];
%VX = (cos(2*pi*(0:q-1)/q-pi)+1).^k;
%VY = (cos(2*pi*(0:p-1)/p-pi)+1).^k;
%VZ = (cos(2*pi*(0:r-1)/r-pi)+1).^k;

[VX,VY,~] = meshgrid(VX,VY,VZ);


D = @(U)D_Forward(U,VX,VY);
Dt = @(Uc)D_Adjoint(Uc,VX,VY,p,q,r);


% high order finite differences
function [Uc] = D_Forward(U,VX,VY)

   
    [a,b,c] = size(U);
    X = fft(U,b,2);
    Y = fft(U,a,1);
    %Z = fft(U,c,3);
    
    X = X.*VX;
    Y = Y.*VY;
    %Z = Z.*VZ;

    X = ifft(X,b,2);
    Y = ifft(Y,a,1);
    %Z = ifft(Z,c,3);
    
    Uc = 2^(1-k)*[X(:),Y(:)];
    



%transpose FD
function dtxy = D_Adjoint(Uc,VX,VY,a,b,c)
    Uc = reshape(Uc,a*b*c,2);
    X = reshape(Uc(:,1),a,b,c);
    Y = reshape(Uc(:,2),a,b,c);
    %Z = reshape(Uc(:,3),a,b,c);
    
    X = fft(X,b,2);
    Y = fft(Y,a,1);
    %Z = fft(Z,c,3);
    
    
    X = X.*conj(VX);
    Y = Y.*conj(VY);
    %Z = Z.*conj(VZ);

    X = ifft(X,b,2);
    Y = ifft(Y,a,1);
    %Z = ifft(Z,c,3);
    
    dtxy = 2^(1-k)*(X(:) + Y(:));
    
    
    
    