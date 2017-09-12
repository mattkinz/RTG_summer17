function [mstack,c] = stack_matrix(stack)

    [m,n,k]=size(stack);
    s = sum(stack);
    s = reshape(s,n,k)';
    [~,c]=min(s);
    c=c';

    y = zeros(m*k*n,3);
    for i =1:n
        for j = 1:k
            y(m*((i-1)*k+j-1)+1:m*((i-1)*k+j),1)=k*(i-1)+j;
            y(m*((i-1)*k+j-1)+1:m*((i-1)*k+j),2)=((k*(i-1)+j-1)*m+1):(k*(i-1)+j)*m; 
            y(m*((i-1)*k+j-1)+1:m*((i-1)*k+j),3)=stack(:,i,j);
        end
    end
    mstack = sparse(y(:,1),y(:,2),y(:,3),n*k,m*n*k);
    
end

