function z = my_cconv3D(x,y)

% my convolution function 

% Written by Toby Sanders @ASU
% Dept. of Statistical & Mathematical Sciences
% 09/22/2016

[m1,n1,k1] = size(x);
[m2,n2,k2] = size(y);

m = max(m1,m2); n = max(n1,n2); k = max(k1,k2);
z = ifftn(fftn(x,[m,n,k]).*conj(fftn(y,[m,n,k])));
