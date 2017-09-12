function [Y] = shearlet_transform_complex(X,shearSys,mode,theta)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016



%compute shearlet coefficients 
switch mode
    case 1
        Y = zeros(size(shearSys.shearlets));
        Xfreq = fftshift(fft2(X.*theta));
        % shearlets coefficients are theoretically computed with a
        % convolution, due to the shifting of the elements
        % Computationally, we can just compute a product of the DFT's
        for i = 1:shearSys.nShearlets
            Y(:,:,i) = ifft2(ifftshift(Xfreq.*conj(shearSys.shearlets(:,:,i))))/shearSys.RMS(i);
        end
        Y = Y*mean(shearSys.RMS);
    case 2
        Y = zeros(size(X,1),size(X,2));
        % transform coefficients into frequency domain
        % similar to case 1
        X = fft2(X);
        for i = 1:shearSys.nShearlets
            Y = Y + fftshift(X(:,:,i)).*shearSys.shearlets(:,:,i)/shearSys.RMS(i);
        end
        Y = ifft2(ifftshift(Y)).*conj(theta);
        Y = Y(:)*mean(shearSys.RMS);
        
end