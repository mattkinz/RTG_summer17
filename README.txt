The goal of these pieces of code is to build a function defined over a square out of simple functions for which we
can easily perform the exact fourier transform to find the fourier coefficients. We can then test the ability of the solvers in irt and 
shared_code_v2.2 to recover the original image (function). We can adjust the pattern over which we "sample" the fourier data and we 
can also add noise to the fourier data. 


You need to run the setup for both irt (setup.m) and the shared_code_v2.2 (getting_started.m).

The code DFT_2D_nonuniform_toby.m is the driver for the whole process. It will call all the other necessary code to perform the experiments.


