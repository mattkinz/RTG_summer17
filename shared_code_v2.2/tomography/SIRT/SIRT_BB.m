function [x,out] = SIRT_BB(stack,angles,recsize,iter,minc,maxc)


%SIRTDEN, SIRTAUTO%


%DESCRIPTION:
    %these function run the SIRT algorithm for tomographic reconstruction.
    
%NOTATION:
    %X = SIRT(stack,angles,recsize,iterations,minc);
    %X = sirtauto(stack,angles,recsize,iterations,minc);

%INPUTS:    
    %stack - the tilt series, where it is assumed the tilt axis is horizontal 
    %and located at the middle of the stack
    %angles - a vector holding the projection angles of the stack, in order,
    %in degrees
    %recsize - the dimension of the reconstruction
    %iter - the number of SIRT iterations
    %minc - the density constraint forcing the reconstruction above that value

    %default values are used for recsize, iter, and minc if they are not
    %specified, therefore one may simply input "sirtden(stack,angles)."
    %the recsize will be set to the detector count, i.e. size(stack,1)
    %50 iterations is default, and no density constraint is used if it is not
    %specified

%OUTPUT: 
    %X - the reconstruction from the input tilt series and other input 
        %parameters

%What is the difference between SIRTDEN and SIRTAUTO?
%It is useful to perform the back projection process over a chunk of slices
%instead of each individual slice to increase stability in the
%reconstructions.  Therefore the "sirtden" function performs the
%backprojection process over the entire input stack.  However, for more
%accuracy in the reconstruction, it is better to do this over smaller
%chunks, like 50 slices instead of 1000, and then put all of the
%chunk reconstructions together in a nice way.  "sirtauto" does this, so
%sirtauto is convenient whenever one is ready to perform the entire SIRT
%reconstruction of all the slices.

%Written by: Toby Sanders @ASU & @PNNL & @USC
%7/11/2014
% Last update: 05/18/2016


% Assign default values to unspecified parameters
if (nargin < 3 || isempty(recsize)), recsize = 512; end
if (nargin < 4 || isempty(iter)), iter = 50; end
if (nargin < 5 || isempty(minc)), constraint = false; minc=[];
else constraint = true; end
if (nargin < 6 || isempty(maxc)), maxconstraint = false; maxc=[];
else maxconstraint = true; end



[W,d1,d2,bb]=getSIRT(angles,recsize,stack);


slices=size(bb,2);
x = zeros(size(W,2),slices);
bbnorm = norm(bb,'fro');
if max(size(recsize))>1
    recsize1 = recsize(1);
    recsize2 = recsize(2);
else
    recsize1 = recsize;
    recsize2 = recsize;
end


errors = zeros(iter,1);
for i = 1:iter
    
    
    
    %compute current difference error
    a1 = bb - W*x;
    errors(i)=norm(a1,'fro')/bbnorm;

    % gradient decent!!
    % d1 and d2 are for normalization
    %g = (W'*(d1.*a1)).*d2; 
    g = W'*a1;
    
    % Compute step length, alpha
    %find constant which minimizes ell_2 error
    if i==1
        a2 = W*g;
        sd = a1'*a2;
        sn = norm(a2,'fro')^2;
        alpha=sd/max(sn,eps);
    else
        xxp = x - xp;
        sd = xxp'*xxp;
        sn = xxp'*(g-gp);
        alpha = abs(sd/sn);
    end
    
    xp = x;gp = g;
    %update
    x = x + g*alpha;
    
    %apply density constraints
    if constraint
        x = max(x,minc);
    end
    if maxconstraint
        x = min(x,maxc);
    end
    
    fprintf('iteration=%i, ||Wx-b||/||b||=%g , ||x||=%g\n',i,errors(i),norm(x,'fro'));
    
    
end

if constraint
    x = max(x,minc);
end
if maxconstraint
    x = min(x,maxc);
end
out.error = errors;
x = reshape(x,recsize1,recsize2,slices);

