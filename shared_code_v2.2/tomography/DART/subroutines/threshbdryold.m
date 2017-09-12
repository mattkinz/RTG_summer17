function [U,set,fixed] = threshbdry(U,threshes,grays)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


num=max(size(threshes));       
[m,n,k]=size(U);
U = reshape(U,m*n,k);
set=cell(k,1);
fixed=cell(k,1);

for i = 1:k
    [s,s2]=sort(U(:,i));
    indexold=0;
    for p = 1:num
        indexup=m*n;
        indexdown=indexold+1;
        index=round((indexup+indexdown)/2);
        while indexup-indexdown~=1
            if s(index)>threshes(p)
                indexup=index;
                index=round((index+indexdown)/2);
            else
                indexdown=index;
                index=round((index+indexup)/2);
            end
        end
        U(s2(indexold+1:index),i)=grays(p);
        indexold=index;
        if indexold==m*n,break;end;
    end
    U(s2(index+1:end),i)=grays(end);
end

U = reshape(U,m,n,k);
mask=zeros(m,n,k,'single');
f = find([diff(U,1,2),zeros(m,1,k)]);
mask(f)=1;
mask(f+m)=1;
f = find([diff(U,1,1);zeros(1,n,k)]);
mask(f)=1;
mask(f+1)=1;
f = find(cat(3,diff(U,1,3),zeros(m,n)));
mask(f)=1;
mask(f+m*n)=1;


mask = reshape(mask,m*n,k);

for i = 1:k
    set{i}=uint32(find(mask(:,i)));
    fixed{i}=uint32(find(mask(:,i)-1));
end







end