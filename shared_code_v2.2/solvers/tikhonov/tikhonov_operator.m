function y = tikhonov_operator(A,D,Dt,dim1,dim2,dim3,datadim,u,mode,lambda)


%A is a function handle that represents the forward operator
%D is the finite difference transform
%Dt is the transpose of D
%dim variables are dimension of the signal
%datadim is the number of data points
%u is the input signal 
%mode is 1 or 2
%lambda is the penalty constant

% Written by Toby Sanders @ASU
% Department of Statistical & Mathematical Sciences
% 05/17/2016

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
        y = y1+lambda*y2;
        
end
        
        