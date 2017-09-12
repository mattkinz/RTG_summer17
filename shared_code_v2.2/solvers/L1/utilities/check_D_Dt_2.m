function [flg,x,y] = check_D_Dt_2(A,N)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


% Check if function handle A has properly defined adjoint
% If it is the true adjoint, then we should have x=y
% flg returns true if Dt is true adjoint and false otherwise
u = rand(N) + 1i*rand(N);
Au = A(u,1);

v = rand(size(Au)) + 1i*rand(size(Au));
Atv = A(v,2);

x = sum(Au(:).*conj(v(:)));
y = sum(u(:).*conj(Atv(:)));

tol = 1e-7;
flg = abs(x-y)<tol;

