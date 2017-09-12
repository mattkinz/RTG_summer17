function [Y] = shearlet_transform(X,shearSys,mode)


% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016



I = find(shearSys.shearletIdxs(:,2)>3);
I = I';
%compute shearlet coefficients 
switch mode
    case 1
        Y = zeros(size(X,1),size(X,2),numel(I));%zeros(size(shearSys.shearlets));
        Xfreq = fftshift(fft2(X));
        for i = I%1:shearSys.nShearlets
            Y(:,:,i) = ifft2(ifftshift(Xfreq.*conj(shearSys.shearlets(:,:,i))))*shearSys.RMS(i);
        end
        Y = Y/mean(shearSys.RMS(I));
    case 2
        Y = zeros(size(X,1),size(X,2));
        %transform coefficients into frequency domain
        X = fft2(X);
        for i = I%1:shearSys.nShearlets
            Y = Y + fftshift(X(:,:,i)).*shearSys.shearlets(:,:,i)*shearSys.RMS(i);
        end
        Y = ifft2(ifftshift(Y));
        Y = Y(:)/mean(shearSys.RMS(I));
        
end