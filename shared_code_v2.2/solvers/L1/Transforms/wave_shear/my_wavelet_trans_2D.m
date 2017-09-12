function [D,Dt] = my_wavelet_trans_2D(wname,levels,p,q)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/22/2016

% these are designed for 2-D problems
[Lo_d,Hi_d] = wfilters(wname);
lo_x = fft(Lo_d,q);
lo_y = fft(Lo_d,p);
hi_x = fft(Hi_d,q);
hi_y = fft(Hi_d,p);

[LOX,LOY] = meshgrid(lo_x,lo_y);
[HIX,HIY] = meshgrid(hi_x,hi_y);




D = @(x)wave_forward(x,LOX,LOY,HIX,HIY,levels,p,q);
Dt = @(c)wave_adjoint(c,LOX,LOY,HIX,HIY,levels,p,q);


function c = wave_forward(U,LOX,LOY,HIX,HIY,levels,p,q)

U = reshape(U,p,q);
Ux = fft(U,q,2);
Uy = fft(U,p,1);

cx = zeros(p,q,levels);
cy = zeros(p,q,levels);
for i = 1:levels
    cx(:,:,i) = ifft(Ux.*HIX,q,2);
    Ux = Ux.*LOX;
    
    cy(:,:,i) = ifft(Uy.*HIY,p,1);
    Uy = Uy.*LOY;
end

c = [cx(:),cy(:)];


function U = wave_adjoint(c,LOX,LOY,HIX,HIY,levels,p,q)

cx = fft(reshape(c(:,1),p,q,levels),q,2);
cy = fft(reshape(c(:,2),p,q,levels),p,1);

Ux = zeros(p,q,levels);
Uy = zeros(p,q,levels);
for i = 1:levels
    Ux(:,:,i) = ifft(cx(:,:,i).*conj(HIX),q,2);
    Uy(:,:,i) = ifft(cy(:,:,i).*conj(HIY),p,1);
    for j = i+1:levels
        cx(:,:,j) = cx(:,:,j).*conj(LOX);
        cy(:,:,j) = cy(:,:,j).*conj(LOY);
    end
end
U = col(sum(Ux,3) + sum(Uy,3));



