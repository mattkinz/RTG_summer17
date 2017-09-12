function x = ARTangle_store(W1,W2,bb,iter,init,numray)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

    N = size(W1,1)/numray;
    
    [m,n]=size(W1);
    if init==0
        x = zeros(1,n,iter+1);
    else
        x = reshape(init,1,n,iter+1);
    end
    p = randi([1 N],iter,1);
    v = zeros(N,2);
    z = ceil(numray/sqrt(n));
    for i = 1:N 
        v(i,1)=(i-1)*numray+1;
        v(i,2) = i*numray-z;
    end
    z = ceil(numray/sqrt(n));
    for i = 1:iter
       list = v(p(i),1)+mod(i,2):z:v(p(i),2)+mod(i,2);
        x(1,:,i+1) = max(0,x(1,:,i)+ (bb(list)-W1(list,:)*x(1,:,i)')'*W2(list,:));
    end

end