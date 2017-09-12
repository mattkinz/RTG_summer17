function stack = weighting2(stack,accuracy)

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


if ~exist('accuracy','var')
    accuracy = 'medium';
elseif ~ischar(accuracy)
    accuracy = 'medium';
end

if  strcmp(accuracy,'good')
    num = 800;
    delta = .025;
elseif strcmp(accuracy,'super')
    num = 2000;
    delta = .01;
elseif strcmp(accuracy,'unbelievable')
    num = 20000;
    delta = .001;
else
    num = 100;
    delta = .1;
end
    


[m,n,k]=size(stack);

[alpha,beta]=min(sum(sum(stack)));
wg = zeros(m,n);

v = 1:k;
v(beta)=0;
r = 1:m;
r = 2/(m-1)*(r-1)-1;
r = cos(pi*r.^2)/2+1/2;
s = zeros(1,k);
        for p = 1:num/10
            rp = r.^(p*delta*10);
            for x =1:n
                    wg(:,x)=rp;
            end
            for i =1:k
                if v(i)
                    if sum(sum(stack(:,:,i).*wg))<alpha
                        v(i)=0;
                        s(i)=(p-1)*10;
                    end
                end
            end
            if sum(v)==0,break;end;
        end
        for i = 1:k
            if v(i),s(i)=(p-1)*10;end
        end
v = 1:k;
v(beta)=0;
    for j = 1:k
        if j~=beta
            for p = 1:10
                rp = r.^((p+s(j))*delta);
                for x =1:n
                        wg(:,x)=rp;
                end
                if sum(sum(stack(:,:,i).*wg))<alpha
                    s(j)=p+s(j);v(i)=0;break;
                end
            end
        end
    end
    for i = 1:k
            if v(i),s(i)=s(i)+10;end
    end
        
    %s=[s(1),s,s(end)];
    %s = imfilter(s,[1/4 1/2 1/4]);
        for i =1:k
                for x = 1:n
                    wg(:,x)=r.^(s(i)*delta);
                end
            stack(:,:,i)=stack(:,:,i).*wg;
        end

        
end
