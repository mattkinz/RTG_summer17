function Fhat = fhat_upper_box_piecewise(kx,ky,A,B1,B2,C1,C2)

%returns the fourier coefficients evaluated over the left-hand side rectangular area.
%We can define the functions 
%f(x,y) = A*squarepulse
%f(x,y) = B1*x
%f(x,y) = B2*y
%f(x,y) = C1*x^2
%f(x,y) = C2*y^2
%
%input
%kx,ky: freqency samples (VECTORS)
%
%output
%Fhat: fourier coefficients
%

z = numel(kx);

%set the boundaries of the rectangle
x1 = -2.5;
x2 = 2.5;
y1 = .5;
y2 = 2.5;

%set the parameters of each function defined on the rectangle
%A = -y1^2;	%f(x,y) = A*squarepulse
%B1 = 0; %f(x,y) = B1*x
%B2 = -2*y1; %f(x,y) = B2*y
%C1 = 0; %f(x,y) = C1*x^2
%C2 = 1; %f(x,y) = C2*y^2

Fhat_square = zeros(z,1);
Fhat_x = zeros(z,1);
Fhat_y = zeros(z,1);
Fhat_xsquared = zeros(z,1);
Fhat_ysquared = zeros(z,1);

%-------------------------------------------------------
%sqaure pulse
if A ~= 0
	for j = 1:z
    	if ky(j) ~= 0 && kx(j) ~= 0
        	P1 = -1/(kx(j)*ky(j)); 
			P2 = exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1);
			P3 = exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1);
			Fhat_square(j) = P1*P2*P3;
    	elseif ky(j) == 0 && kx(j) ~= 0
        	P1 = 1i/kx(j); 
			P2 = exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1);
			P3 = y2-y1;
			Fhat_square(j) = P1*P2*P3;
		elseif ky(j) ~= 0 && kx(j) == 0
    	   	P1 = 1i/ky(j); 
		   	P2 = x2-x1;
	    	P3 = exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1);
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
	for j = 1:z
    	if ky(j) ~= 0 && kx(j) ~= 0
        	P1 = 1i/ky(j) * (exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1));
			P2 = (1i*x2/kx(j) + 1/kx(j)^2) * exp(-1i*kx(j)*x2); 
			P3 = (1i*x1/kx(j) + 1/kx(j)^2) * exp(-1i*kx(j)*x1); 
			Fhat_x(j) = P1*(P2-P3);
	    elseif ky(j) == 0 && kx(j) ~= 0
    	    P1 = y2 - y1;
			P2 = (1i*x2/kx(j) + 1/kx(j)^2) * exp(-1i*kx(j)*x2); 
			P3 = (1i*x1/kx(j) + 1/kx(j)^2) * exp(-1i*kx(j)*x1); 
			Fhat_x(j) = P1*(P2-P3);
	    elseif ky(j) ~= 0 && kx(j) == 0
        	P1 = 1i/ky(j) * (exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1));
			P2 = 1/2 * (x2^2 - x1^2);
			Fhat_x(j) =P1*P2;
	    else 
    	    P1 = y2 - y1;
			P2 = 1/2 * (x2^2 - x1^2);
			Fhat_x(j) =P1*P2;
	    end
	end
	Fhat_x = B1*Fhat_x;
end
%--------------------------------------------------------
%f(x,y) = B2*y
if B2 ~= 0
	for j = 1:z
    	if kx(j) ~= 0 && ky(j) ~= 0
        	P1 = 1i/kx(j) * (exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1));
			P2 = (1i*y2/ky(j) + 1/ky(j)^2) * exp(-1i*ky(j)*y2); 
			P3 = (1i*y1/ky(j) + 1/ky(j)^2) * exp(-1i*ky(j)*y1); 
			Fhat_y(j) = P1*(P2-P3);
	    elseif kx(j) == 0 && ky(j) ~= 0
    	    P1 = x2 - x1;
			P2 = (1i*y2/ky(j) + 1/ky(j)^2) * exp(-1i*ky(j)*y2); 
			P3 = (1i*y1/ky(j) + 1/ky(j)^2) * exp(-1i*ky(j)*y1); 
			Fhat_y(j) = P1*(P2-P3);
	    elseif kx(j) ~= 0 && ky(j) == 0
        	P1 = 1i/kx(j) * (exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1));
			P2 = 1/2 * (y2^2 - y1^2);
			Fhat_y(j) =P1*P2;
	    else 
    	    P1 = x2 - x1;
			P2 = 1/2 * (y2^2 - y1^2);
			Fhat_y(j) =P1*P2;
	    end
	end
	Fhat_y = B2*Fhat_y;
end
%-------------------------------------------------------------
%f(x,y) = C1*x^2
if C1 ~= 0
	for j = 1:z
    	if ky(j) ~= 0 && kx(j) ~= 0
        	P1 = 1i/ky(j) * (exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1));
			P2 = (1i*x2^2/kx(j) + 2*x2/kx(j)^2 + 2/(1i*kx(j)^3)) * exp(-1i*kx(j)*x2); 
			P3 = (1i*x1^2/kx(j) + 2*x1/kx(j)^2 + 2/(1i*kx(j)^3)) * exp(-1i*kx(j)*x1); 
			Fhat_xsquared(j) = P1*(P2-P3);
	    elseif ky(j) == 0 && kx(j) ~= 0
    	    P1 = y2 - y1;
			P2 = (1i*x2^2/kx(j) + 2*x2/kx(j)^2 + 2/(1i*kx(j)^3)) * exp(-1i*kx(j)*x2); 
			P3 = (1i*x1^2/kx(j) + 2*x1/kx(j)^2 + 2/(1i*kx(j)^3)) * exp(-1i*kx(j)*x1); 
			Fhat_xsquared(j) = P1*(P2-P3);
	    elseif ky(j) ~= 0 && kx(j) == 0
        	P1 = 1i/ky(j) * (exp(-1i*ky(j)*y2) - exp(-1i*ky(j)*y1));
			P2 = 1/3 * (x2^3 - x1^3);
			Fhat_xsquared(j) =P1*P2;
	    else 
    	    P1 = y2 - y1;
			P2 = 1/3 * (x2^3 - x1^3);
			Fhat_xsquared(j) =P1*P2;
	    end
	end
	Fhat_xsquared = C1 * Fhat_xsquared;
end
%-------------------------------------------------------------------
%f(x,y) = C2*y^2
if C2 ~= 0
	for j = 1:z
    	if kx(j) ~= 0 && ky(j) ~= 0
        	P1 = 1i/kx(j) * (exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1));
			P2 = (1i*y2^2/ky(j) + 2*y2/ky(j)^2 + 2/(1i*ky(j)^3)) * exp(-1i*ky(j)*y2); 
			P3 = (1i*y1^2/ky(j) + 2*y1/ky(j)^2 + 2/(1i*ky(j)^3)) * exp(-1i*ky(j)*y1); 
			Fhat_ysquared(j) = P1*(P2-P3);
	    elseif kx(j) == 0 && ky(j) ~= 0
    	    P1 = x2 - x1;
			P2 = (1i*y2^2/ky(j) + 2*y2/ky(j)^2 + 2/(1i*ky(j)^3)) * exp(-1i*ky(j)*y2); 
			P3 = (1i*y1^2/ky(j) + 2*y1/ky(j)^2 + 2/(1i*ky(j)^3)) * exp(-1i*ky(j)*y1); 
			Fhat_ysquared(j) = P1*(P2-P3);
	    elseif kx(j) ~= 0 && ky(j) == 0
        	P1 = 1i/kx(j) * (exp(-1i*kx(j)*x2) - exp(-1i*kx(j)*x1));
			P2 = 1/3 * (y2^3 - y1^3);
			Fhat_ysquared(j) =P1*P2;
	    else 
    	    P1 = x2 - x1;
			P2 = 1/3 * (y2^3 - y1^3);
			Fhat_ysquared(j) =P1*P2;
	    end
	end
	Fhat_ysquared = C2 * Fhat_ysquared;
end
%-----------------------------------------------------------------------

Fhat = (Fhat_square + Fhat_x + Fhat_y + Fhat_xsquared + Fhat_ysquared) / (2*pi)^2;
Fhat = reshape(Fhat,sqrt(numel(ky)),sqrt(z));
end
