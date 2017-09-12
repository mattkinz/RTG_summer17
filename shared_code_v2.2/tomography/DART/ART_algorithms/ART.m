function x = ART(W1,W2,bb,iter,init)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

    [m,n]=size(W1);
    if init==0
        x = zeros(1,n);
    else
        x = reshape(init,1,n);
    end
    p = randi([1 m],iter,1);
    max(p)
    for i = 1:iter
        x = x+ (bb(p(i))-W1(p(i),:)*(x'))*W2(p(i),:);
    end

end