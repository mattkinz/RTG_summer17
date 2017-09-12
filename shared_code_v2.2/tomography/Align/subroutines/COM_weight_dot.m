function f = COM_weight_dot(stack,w);

[m,n,k]=size(stack);
f = zeros(n,k);
for i = 1:n
    for j = 1:k
        f(i)=stack(:,i,j)'*w((i+j-2)*m+1:(i+j-1)*m);
    end
end

end