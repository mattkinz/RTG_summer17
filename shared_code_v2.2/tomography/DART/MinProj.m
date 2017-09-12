
function [p] = MinProj(stack,angles,U,thresh)


%This function computes the gray values which minimize the projection
%error from the given input data.  These minimizing gray values are
%found individually for slices of the reconstruction.  All of the
%information is stored in "p" and described below.  If it is not clear,
%this function additionally prints an output for clarity.
    
%INPUTS:
    %angles - an order list of the projections angles
    %or angles may simply be set to the projection matrix
    %U - the original sirt or tvm, or otherwise reconstruction
    %stack - the aligned tilt series, where it is assumed the tilt axis is
        %horizontal and in the middle of the stack
    %threshes - set of thresholds (usually one or two), estimated for
        %the input reconstruction
%Output:
    %p-     
    %In p(j,1) is the slice number corresponding to the data in p(j,2:end).
    %The minimum projection gray values for slice p(j,1) will be stored in 
    %p(:,2:numgrays+1) where numgrays is the number of gray values, which 
    %is one greater than the length of the input vector thresh.    
    %In p(j,numgrays+2) is the projection error for slice p(j,1) before 
    %segmentation and in p(j,numgrays+3) is the projection error 
    %after segmentation determined by the minimum projection gray values.
    
    %The second set of gray values does not assume the lowest gray level to
    %be zero, and is hence the lowest gray level is also output
    %The error from using the gray values is in p(end), and obviously
    %p(end)<=p(1).  It is suggested to use this second set of gray levels


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



increment=10;

if min(size(angles))==1
    W = radonmatrix(angles,size(U,1),size(stack,1),0);
else
    W=angles;
end
grays = max(size(thresh));

[m,n,k]=size(U);
p=zeros(ceil(k/increment),grays+4);
if k~=size(stack,2)
    error('stack and reconstruction do not match');
end


for j = 1:ceil(k/increment)
    
Ut = U(:,:,(j-1)*increment+1);
stackt = stack(:,(j-1)*increment+1,:);



if size(stackt,3)~=1
    stackt = reshape(stackt,...
        size(stackt,1)*size(stackt,2)*size(stackt,3),1);
end

Ut = reshape(Ut,m*n,1);
[s,s2]= sort(Ut);
index=1;
c=1;
counters = zeros(grays+1,1);
for i = 1:m*n
    if s(i)>=thresh(c)
	counters(c)=i-index;
	c = c+1;
	index=i;
	if c == grays+1;
	    counters(c)=m*n-index+1;
	    break;
    end
    end
end


P = zeros(size(W,1),grays+1);
for i = 1:grays+1
  P(:,i)=sum(W(:,s2(sum(counters(1:i-1))+1:sum(counters(1:i)))),2);
end


%p(j,2:grays+1) = pinv(P(:,2:end))*stackt;
p(j,2:grays+2)=pinv(P)*stackt;
p(j,1)=(j-1)*increment+1;
nn=norm(stackt,'fro');
p(j,grays+4)=norm(P*p(j,2:grays+2)'-stackt,'fro')/nn;
p(j,grays+3)=norm(W*Ut-stackt,'fro')/nn;
if j==1
fprintf('____slice_____')

for i = 1:grays+1
fprintf('gray%d_____',i);
end
fprintf('error1____');
fprintf('error2\n');
end

disp(p(j,:))

end

end

