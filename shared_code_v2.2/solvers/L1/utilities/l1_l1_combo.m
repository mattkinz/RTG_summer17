function y = l1_l1_combo(A,D,Dt,mu,dims1,dims2,dims3,u,mode)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

switch mode
    case 1
        y = mu*A(u,1);
        y = [y;col(D(reshape(u,dims1)))];
    case 2
        %whos
        y = mu*A(u(1:dims3),2);
        y = y + col(Dt(reshape(u(dims3+1:end),dims2)));
end