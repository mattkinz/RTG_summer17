function [set,fixed,region] = regionfind(U,Uchunk,bb,W,thresh,tol)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


[m,n,k]=size(U);

mask = zeros(m,n,k,'single');

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
Uchunk = reshape(Uchunk,m*n,k);
U = reshape(U,m*n,k);
region =nnz(mask);
forwardprojdiff = abs(W*(Uchunk-U));
for i = 1:k
    S = max((forwardprojdiff(:,i))./(bb(:,i)+10)-1/2,0);
    y = sum(W(logical(S),:));
    T = find(max(y-1,0));
    for l = 1:size(T,2)
        for j = 1:size(thresh)           
            if thresh(j)-tol<Uchunk(T(l),i) && Uchunk(T(l),i)<thresh(j)+tol
                mask(T(l),i)=1;   
                break;
            end
        end
    end
    %mask(logical(max(y-30,0)),i)=1;
end
region = nnz(mask)-region
set = cell(k,1);
fixed = cell(k,1);
for i = 1:k
    set{i}=uint32(find(mask(:,i)));
    fixed{i}=uint32(find(mask(:,i)-1));
end


end


