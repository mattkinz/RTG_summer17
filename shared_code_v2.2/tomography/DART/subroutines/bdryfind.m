function [set,fixed] = bdryfind(U,thresh,grays)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



[m,n,k]=size(U);
dx = diff(U,1,1);
dy = diff(U,1,2);
dz = diff(U,1,3);


mask=zeros(m,n,k,'single');
for j = 1:m-1
    for i=1:k;for l = 1:n
            if dx(j,l,i)
                mask(j,l,i)=1;mask(j+1,l,i)=1;
            end
        end;end
end

for l = 1:n-1
    for i =1:k;for j=1:m
            if dy(j,l,i)
                mask(j,l,i)=1;mask(j,l+1,i)=1;
            end
        end;end;
end
               
for i = 1:k-1
    for j = 1:m; for l = 1:n
            if dz(j,l,i)
              mask(j,l,i)=1;mask(j,l,i+1)=1;
            end
    end;end;
end


mask = reshape(mask,m*n,k);
U = reshape(U,m*n,k);
set = cell(k,1);
fixed = cell(k,1);

for i = 1:k
    set{i}=uint32(find(mask(:,i)));
    fixed{i}=uint32(find(mask(:,i)-1));
end





end