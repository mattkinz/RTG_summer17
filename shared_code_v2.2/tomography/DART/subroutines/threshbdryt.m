function [U,set,fixed] = threshbdry(U,threshes,grays)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


num=max(size(threshes));       
[m,n,k]=size(U);
U = reshape(U,m*n,k);
if num>1
    a = 9/8*threshes(1)-1/8*threshes(end);
    b = (threshes(end)-a)/.9;
    threshes = (threshes-a)/b;
    U = (U-a)/b;
else
    U = U-threshes+1/2;
    threshes=1/2;
end
set=cell(k,1);
fixed=cell(k,1);
Ut = cell(num,1);

for i = 1:num
    Ut{i}=im2bw(U,threshes(i));
end
U = zeros(m*n,k)+grays(1);
for i =1:num
    U = U + double(Ut{i})*(grays(i+1)-grays(i));
end

U = reshape(U,m,n,k);
mask=zeros(m,n,k,'single');
f = find([diff(U,1,2),zeros(m,1,k)]);
mask(f)=1;
mask(f+m)=1;
mask(f+m-1)=1;mask(min(f+m+1,m*n*k))=1;
f = find([diff(U,1,1);zeros(1,n,k)]);
mask(f)=1;
mask(f+1)=1;
mask(min(f+1+m,m*n))=1;
mask(max(1,f+1-m))=1;
f = find(cat(3,diff(U,1,3),zeros(m,n)));
mask(f)=1;
mask(f+m*n)=1;



mask = reshape(mask,m*n,k);

for i = 1:k
    set{i}=uint32(find(mask(:,i)));
    fixed{i}=uint32(find(mask(:,i)-1));
end







end