function [U,set,fixed] = thresh_tol(U,grays,tol)


%The tolerance is around the gray values, not the thresholds
[m,n,k]=size(U);
U = reshape(U,m*n,k);
num = max(size(grays))-1;
r = zeros(2*num,1);
r(1)=grays(1)+tol;
r(end)=grays(end)-tol;
fixed = cell(k,1);
set = cell(k,1);
for i = 1:num-1
    r(2*i)=grays(i+1)-tol;
    r(2*i+1)=grays(i+1)+tol;
end

for i = 1:k
    [s,s2]=sort(U(:,i));
    c=1;
    indexold=0;
    locater=zeros(2*num+1,1);
    for j = 1:2*num
        indexdown = indexold+1;
        indexup=m*n;
        while ~le(indexup-indexdown,1)
            index = round(mean([indexup indexdown]));
            if s(index)<r(c)
                indexdown=index;
            else
                indexup=index;
            end
        end
        locater(c)=index;
        c = c+1;
        indexold=index;
        if indexold==m*n
            break;
        end
    end
    locater(c:end)=m*n;
    U(s2(1:locater(1)),i)=grays(1);
    U(s2(locater(end-1):end),i)=grays(end);
    for j = 1:num-1
        U(s2(locater(2*j):locater(2*j+1)),i)=grays(j+1);
    end
    set{i}=[];
    fixed{i}=s2(1:locater(1));

    for p = 1:num
        fixed{i}=[fixed{i};s2(locater(2*p):locater(2*p+1))];
        set{i}=[set{i};s2(locater(2*p-1)+1:locater(2*p)-1)];
    end
end
U = reshape(U,m,n,k);