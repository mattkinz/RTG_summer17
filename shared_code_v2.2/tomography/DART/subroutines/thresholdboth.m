
function [U,set,fixed] = thresholdboth(U,threshes,grays,tol)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



num=max(size(threshes));
if max(size(tol))==1
    tol=tol*ones(num,1);
end
for i = 1:num
    if tol(i)>=threshes(i)-grays(i) || tol(i)>=grays(i+1)-threshes(i)
        error('the tolerance specified is too large\n');
    end
end
 
        
[m,n,k]=size(U);
U = reshape(U,m*n,k);

set=cell(k,1);
fixed=cell(k,1);
mask = zeros(m*n,k,'single');
for i = 1:k
    l1 = zeros(num,2);
    l2 = zeros(num+2,1);
    l2(1)=1;l2(end)=m*n;
    [s,s2]=sort(U(:,i));
    c=1;
    check=0;
    for j = 1:m*n
        if check==0
            if threshes(c)-s(j)<tol(c)
                l1(c,1)=j;
                if s(j)<threshes(c)
                    check=1;
                else
                    check=2;
                    c=c+1;
                    l2(c)=j-1;
                end
            end
        elseif check==1
            if s(j)>threshes(c)
                c=c+1;
                l2(c)=j-1;             
                if s(j)-threshes(c-1)<tol
                    check=2;                    
                else
                    check=0;
                    l1(c-1,2)=j-1;
                    if c == num+1, break;end;
                end                   
            end
        elseif check==2
            if s(j)-threshes(c-1)>=tol(c-1)
                l1(c-1,2)=j-1;
                check=0;
                if c==num+1, break;end;
            end
        end        
    end
    for l = 1:num
        if l1(l,1)~=0 && l1(l,2)~=0
            mask(s2(l1(l,1):l1(l,2)),i)=1;
        elseif l1(l,1)~=0
            mask(s2(l1(l,1):m*n),i)=1;
    end
    end

    for l = 1:num+1
        if l2(l+1)~=0
            U(s2(l2(l)+1:l2(l+1)),i)=grays(l);
        else
            U(s2(l2(l)+1:end),i)=grays(l);
            break;
        end
    end
end

U = reshape(U,m,n,k);
mask = reshape(mask,m,n,k);
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
            
            