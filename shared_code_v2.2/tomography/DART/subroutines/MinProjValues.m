
function [p] = MinProjValues(angles,U,stack,threshes)


%angles are the projection angles
%U is the estimated solution
%bb is the data, stacked or not
%threshes should be a cell contained the threshold values to be tested
%for example threshes{2} should be a vector holding testing values for the 2nd threshhold
%all the combinations of thresh values will be paired
%for each pairing, the gray values which minimize the thresholding for that pairing will be found
%the residual for that thresholding is then computed
%we want to choose the thresholding that gives the smallest residual
%p(:,1) holds the residuals
%p(:,2:grays+1) holds the grays values that minimize the thresholding 
%for the thresholding values listed in p(:,grays+2:2*grays+1)
%p(:,2*grays+2:3*grays+2) holds the gray values without assuming the low is zero
%p(:,end) holds the residual for this case (low gray not nec zeros)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



if min(size(angles))==1
    W = radonmatrix(angles,size(U,1),size(stack,1),0);
else
    W=angles;
end
grays = max(size(threshes));
for i = 1:grays
    t(i) = max(size(threshes{i}));
end
p=zeros(prod(t),3*grays+3);

if size(bb,2)~=1 || size(bb,3)~=1
    bb=reshape(bb,size(bb,1)*size(bb,2)*size(bb,3),1);
end
[m,n]=size(U);
U = reshape(U,m*n,1);
[s,s2]= sort(U);
thresh = zeros(1,grays);
for k = 1:prod(t)
    for j = 1:grays
	thresh(j) =threshes{j}(mod(ceil(k/prod(t(j+1:end)))-1,t(j))+1);
    end
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


p(k,2:grays+1) = pinv(P(:,2:end))*bb;
p(k,2*grays+2:3*grays+2)=pinv(P)*bb;
p(k,1)=norm(P(:,2:end)*p(k,2:grays+1)'-bb);
p(k,3*grays+3)=norm(P*p(k,2*grays+2:3*grays+2)'-bb);
p(k,grays+2:2*grays+1)=thresh;
end

p(:,1)=p(:,1)/1000;
p(:,end)=p(:,end)/1000;
fprintf('Original solution differenc: %d\n',norm(W*U-bb));
fprintf('____Resid1____')
for i = 1:grays
fprintf('gray%d____',i);
end
for i = 1:grays
fprintf('thresh%d____',i);
end
for i = 1:grays+1
fprintf('gray%d_____',i);
end
fprintf('Resid2\n');
disp(p)
end

