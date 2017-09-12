function [stack,shifts] = corrx(stack)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



stackw = weighting2(stack);

[m,n,k]=size(stack);
shifts = zeros(k,1);
X = zeros(k,n);
for i =1:k
    X(i,:)=sum(stackw(:,:,i));
end
a = (n+1)/2;
w = zeros(k,n);
for x = 1:n
    w(:,x)=(1-(x-a)^2/2/a^2);
end
X = X.*w;
imagesc(X);figure(1);
Y = mean(X);

coarse=-20:5:20;
startslice = ceil(k/2);
for i =startslice-1:-1:1
    maxshift = corrsinglex(X(i+1,:),X(i,:),coarse);
    finex = maxshift-3:maxshift+3;
    shifts(i) = corrsinglex(X(i+1,:),X(i,:),finex);
end
for i = startslice+1:k
    maxshift = corrsinglex(X(i-1,:),X(i,:),coarse);
    finex = maxshift-3:maxshift+3;
    shifts(i) = corrsinglex(X(i-1,:),X(i,:),finex);
end

for i = 1:startslice-1
    stack(:,:,i)=circshift(stack(:,:,i),[0,sum(shifts(i:startslice))]);
end
for i = startslice+1:k
    stack(:,:,i)=circshift(stack(:,:,i),[0,sum(shifts(startslice:i))]);
end



    

