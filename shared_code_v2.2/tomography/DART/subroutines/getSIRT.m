function [W,d1,d2,bb]=getSIRT(angles,recsize,stack)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014
% Last update: 05/2015


%determine the number of slices and build projection matrix if angles are
%specified instead of matrix
[slices]=size(stack,2);
if size(angles,1)==1 || size(angles,2)==1
    if size(stack,3)==1
        numray = size(stack,1)/max(size(angles));
    elseif size(stack,3)~=max(size(angles))
        error('Number of input angles doesnt match the number of projection images');
    else
        numray = size(stack,1);
    end
    W = radonmatrix(angles,recsize,numray,0);
else
    W = angles;
end


%get normalizations for smearing
d1 = 1./sum(W,2);
d1 = repmat(d1,1,slices);
d2 = 1./sum(W,1)';
d2 = repmat(d2,1,slices);


%reshape the stack for computations
if size(stack,3)~=1
    bb = zeros(size(stack,1)*size(stack,3),slices);
    for i = 1:slices
        bb(:,i)=reshape(stack(:,i,:),size(stack,1)*size(stack,3),1);
    end
else
    bb=stack;
end