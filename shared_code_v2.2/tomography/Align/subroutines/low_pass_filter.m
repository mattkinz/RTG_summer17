function stack = low_pass_filter(stack,level)

%applies low pass fourier filter to a stack of images
%level should be between 0 and 1, telling which percentage of fourier
%values "pass."  For example, if level is .2, then the smallest 20 percent
%of fourier values will pass, and the rest are set to 0.



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


    [m,n,k]=size(stack);
    cut_off = sqrt(m*m/4+n*n/4)*level;
    for i =1:k
        stack(:,:,i)=fftshift(fftn(stack(:,:,i)));
    end
    
    for i = 1:m
        for j =1:n
            dist = sqrt((m/2-i)^2+(n/2-j)^2);
            if dist>cut_off
                stack(i,j,:)=0;
            end
        end
    end
    
    
    for i =1:k
        stack(:,:,i)=real(ifftn(ifftshift(stack(:,:,i))));
    end
    
end
    
