function Fhat = fhat_left_background_piecewise(kx,ky)

%returns the fourier coefficients evaluated over the left-hand side rectangular area
%input
%kx,ky: freqency samples (VECTORS)
%
%output
%Fhat: fourier coefficients
%

z = numel(kx);
x1 = -pi;
x2 = -2.5;
y1 = -pi;
y2 = pi;

A = 1;	%f(x,y) = A*squarepulse
B1 = 1; %f(x,y) = B1*x
B2 = 1; %f(x,y) = B2*y
C1 = 1; %f(x,y) = C1*x^2
C2 = 1; %f(x,y) = C2*y^2

if A ~= 0
%-------------------------------------------------------
%sqaure pulse
	Fhat_square = zeros(z,1);
	for j = 1:z
    	if ky(j) ~= 0 && kx(j) ~= 0
        	P1 = -1/(kx(j)*ky(j)); 
			P2 = exp(-1i*kx(j)*x2) - exp(1i*kx(j)*x1);
			P3 = exp(-1i*ky(j)*y2) - exp(1i*ky(j)*y1);
			Fhat_square(j) = P1*P2*P3;
    	elseif ky(j) == 0 && kx(j) ~= 0
        	P1 = 1i/kx(j); 
			P2 = exp(-1i*kx(j)*x2) - exp(1i*kx(j)*x1);
			P3 = 2;
			Fhat_square(j) = P1*P2*P3;
	    elseif ky(j) ~= 0 && kx(j) == 0
    	    P1 = 1i/ky(j); 
			P2 = 2;
			P3 = exp(-1i*ky(j)*y2) - exp(1i*ky(j)*y1);
			Fhat_square(j) = P1*P2*P3;
	    else 			%both kx(j) and ky(j) are zero
    	    P1 = x2-x1;
			P2 = y2-y1;
			Fhat_square(j) = P1*P2;
    	end
	end
	Fhat_square = A*Fhat_square;
end
%-------------------------------------------------------
%f(x,y) = B1*x
if B1 ~= 0
	Fhat_x = zeros(z,1);
	for j = 1:z
    	if ky(j) ~= 0 && kx(j) ~= 0
        	P1 = 1i/ky(j) * (exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1));
			P2 =
			P3 =
			Fhat_x(j) = ;
	    elseif ky(j) == 0 && kx(j) ~= 0
    	    P1 = 
			P2 = 
			P3 = 
			Fhat_x(j) = ;
	    elseif ky(j) ~= 0 && kx(j) == 0
    	    P1 = 
			P2 =
			P3 =
			Fhat_x(j) =;
	    else 
    	    P1 = 
			P2 =
			P3 =
			Fhat_x(j) = ;
	    end
	end
end




%at the end make sure to divide by (2*pi)^2
