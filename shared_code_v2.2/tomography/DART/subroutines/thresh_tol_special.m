function [U,set,fixed] = thresh_tol_special(U,threshes,grays)


%The tolerance is around the gray values, not the thresholds
[m,n,k]=size(U);
U = reshape(U,m*n,k);
U2 = zeros(size(U));
num = max(size(threshes));
fixed = cell(k,1);
set = cell(k,1);
i1 = zeros(k,1);
i2 = zeros(k,1);

for i = 1:k
    [s,s2]=sort(U(:,i));

    indexold=0;
    indexoldo=0;
    for j = 1:2
        indexdown = indexold+1;
        indexup=m*n;
        while indexup-indexdown~=1
                index = round(mean([indexup indexdown]));
                if s(index)<threshes(j)
                    indexdown=index;
                else
                    indexup=index;
                end
        end    
      
        if j==1
            U(s2(indexoldo+1:index),i)=grays(j);
            i1(i) = index;
        else
            U(s2(index:end),i)=grays(j);
            i2(i)=index;
        end
        indexold=index;
        if indexold==m*n
            U(s2(index:end),i)=grays(j);
            i2=index;
            break;
        end
    end
    U(s2(index+1:m*n),i)=grays(end);
    set{i}=zeros(m*n,1);
    set{i}(s2(i1(i):i2(i)))=s2(i1(i):i2(i));
    fixed{i}=(1:m*n)';
    U2(:,i)=U(:,i);
    U2(s2(i1(i):i2(i)),i)=mean(grays);
end
U = reshape(U,m,n,k);


U2 = reshape(U2,m,n,k);
mask=zeros(m,n,k,'single');
f = find([diff(U2,1,2),zeros(m,1,k)]);
mask(f)=1;
mask(f+m)=1;
f = find([diff(U2,1,1);zeros(1,n,k)]);
mask(f)=1;
mask(f+1)=1;
f = find(cat(3,diff(U2,1,3),zeros(m,n)));
mask(f)=1;
mask(f+m*n)=1;


mask = reshape(mask,m*n,k);

for i = 1:k
    set{i}(find(mask(:,i)))=uint32(find(mask(:,i)));
    set{i}=find(set{i});
    fixed{i}(set{i})='';
end



