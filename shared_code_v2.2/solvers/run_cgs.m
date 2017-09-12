function X = run_cgs(A,b,tol,maxit)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 04/18/2016



    if ~isa(A,'function_handle')
        A = @(u,mode) f_handleA(A,u,mode);
    end

    if ~exist('tol','var')
        tol = [];
    end

    if ~exist('maxit','var')
        maxit = [];
    end


    B = @(x)A(A(x,1),2);
    b = A(b,2);

    X = cgs(B,b,tol,maxit);

end