function F = f_background_top_piecewise(x,y,A 

%returns the function values evaluated over the specified rectangular area.
%We can define the functions 
%f(x,y) = A*squarepulse
%f(x,y) = B1*x
%f(x,y) = B2*y
%f(x,y) = C1*x^2
%f(x,y) = C2*y^2
%
%input
%x,y: freqency samples (VECTORS)
%A,B1,B2,C1,C2: function parameters (SCALARS)
%
%output
%F: function values
%
F = zeros(numel(ky),numel(kx)) 	%y values down the rows and x values across the columns



%Now do each case for the different functions
