function stack = weighting(stack)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


[m,n,k]=size(stack);

[alpha,beta]=min(sum(sum(stack)));
a=(m+1)/2;
b=(n+1)/2;
wg = zeros(m,n);

v = 1:k;
v(beta)=0;
s = zeros(1,k);
        for p = 1:30
            for y =1:m
                    wg(y,:)=sqrt((2-(y-a)^2/a^2)/2)^p;
            end
            for i =1:k
                if v(i)
                    if sum(sum(stack(:,:,i).*wg))<alpha
                        v(i)=0;
                        s(i)=p;
                    end
                end
            end
        end
        for i = 1:k
            if v(i)
                s(i)=p;
            end
        end
        
        ff = [1/4 1/2 1/4];
        s = [s(1),s,s(end)];
        s = imfilter(s,ff);
        for i =1:k
            for y =1:m
                for x = 1:n
                    wg(y,x)=sqrt((2-(y-a)^2/a^2)/2)^s(i+1);
                end
            end
            stack(:,:,i)=stack(:,:,i).*wg;
        end

        
end
