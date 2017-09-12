function [normDU,dU] = diff_norm_k(U,p,k)

%Written by: Toby Sanders
%Comp. & Applied Math Dept., Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016


%Finds the \ell_p norm of the kth order finite difference of U

if round(k)==k
    [D,~] = FD3D(k,size(U,1),size(U,2),size(U,3));
else
    [D,~] = FFD3D(k,size(U,1),size(U,2),size(U,3));
end

dU = D(U);
normDU = norm(dU(:),p);
