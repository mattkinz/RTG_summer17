function [A,bpad,scl,opts] = tikhonov_setup(A,b,dim1,dim2,dim3,opts)

% Written by Toby Sanders @ASU
% Department of Statistical & Mathematical Sciences
% 05/17/2016

[A,bpad,scl,opts] =...
    get_tikhonov_local(A,b,dim1,dim2,dim3,opts);






function [A,bpad,scl,opts] = get_tikhonov_local(A,b,dim1,dim2,dim3,opts)
%order of the finite difference operator
if ~isfield(opts,'order')
    opts.order = 0;
end

%lambda is the penalty weight for the second norm
if ~isfield(opts,'lambda')
    opts.lambda = .15;%/sqrt(factorial(2*(opts.order))/factorial(opts.order)^2);
    fprintf('lambda not specified, set to default\n');
    fprintf('Default lambda = %f\n',opts.lambda);
end


%scale lambda according to the order
if ~isfield(opts,'scale_lam')
    opts.scale_lam = true;
end


%here we use the fact that for a simple jump signal f, from alpha to beta, 
%we have
% ||D_{k+1} f||_2^2 = (alpha-beta)^2 (2k)!/(k!)^2
if opts.scale_lam
    opts.lambda = opts.lambda/...
        sqrt(factorial(2*(opts.order-1))/factorial(opts.order-1)^2);
end

%maximum iterations for iterative method
%if ~isfield(opts,'maxit')
%    opts.maxit=50;
%end
%maxit = opts.maxit;

%tolerance for convergence
%if ~isfield(opts,'tol')
%    opts.tol = 10^(-3);
%end
%tol = opts.tol;

%initial guess
%if ~isfield(opts,'x0')
%    opts.x0 = zeros(dim1,dim2,dim3);
%end
%x0 = opts.x0(:);



%scaling of the operator A and vector b
if isfield(opts,'scale_A')
    if ~islogical(opts.scale_A)
        error('opts.scale_A should be true or false.');
    end
else
    opts.scale_A = true;
end


if isfield(opts,'scale_b')
    if ~islogical(opts.scale_b)
        error('opts.scale_b should be true or false.');
    end
else
    opts.scale_b = true;
end

Amain = A;
if ~isa(A,'function_handle')
    Amain = @(u,mode) f_handleA(Amain,u,mode);
end

% check scaling A
if opts.scale_A
    [Amain,b] = ScaleA(dim1*dim2*dim3,Amain,b);
end

% check scaling b
scl = 1;
if opts.scale_b
    [b,scl] = Scaleb(b);
end
datadim = size(b,1);


%pad b for the regularization term
bpad = [b;zeros(3*dim1*dim2*dim3,1)];

%define the finite difference operators
[D,Dt] = FD3D(opts.order);

%build the operator [A;D] and [A^t + Dt]
A = @(u,mode)tikhonov_operator_local(Amain,D,Dt,dim1,dim2,dim3,datadim,u,mode,opts.lambda);





function y = tikhonov_operator_local(A,D,Dt,dim1,dim2,dim3,datadim,u,mode,lambda)


%A is a function handle that represents the forward operator
%D is the finite difference transform
%Dt is the transpose of D
%dim variables are dimension of the signal
%datadim is the number of data points
%u is the input signal 
%mode is 1 or 2
%lambda is the penalty constant
switch mode
    case 1
        y1 = A(u(:),1);
        [x,y,z] = D(reshape(u,dim1,dim2,dim3));
        y2 = lambda*[x(:);y(:);z(:)];
        y=[y1;y2];
        
    case 2
        y1 = A(u(1:datadim),2);
        dd = dim1*dim2*dim3;
        x0 = reshape(u(datadim+1:datadim+dd),dim1,dim2,dim3);
        y0 = reshape(u(datadim+dd+1:datadim+2*dd),dim1,dim2,dim3);
        z0 = reshape(u(datadim+2*dd+1:end),dim1,dim2,dim3);
        y2 = Dt(x0,y0,z0);
        y = [y1+lambda*y2];
        
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        