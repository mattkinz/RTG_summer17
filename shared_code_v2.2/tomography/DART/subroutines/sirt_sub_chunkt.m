function x = sirt_sub_chunk(W,Wn,bb,iter,x,set,c,minc,grays)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



slices = max(size(set));
xnew = cell(slices,1);

for i = 1:slices
    subdim=max(size(set{i}));
    xnew{i}=zeros(subdim,1);
end

for j = 1:slices 
    Wj = W(:,set{j});
for i = 1:iter
    %sn=0;sd=0;
          
        a1=bb(:,j)-Wj*x{j};
        xnew{j}=(a1'*Wn(:,set{j}))'.*c{j};
        a2=W(:,set{j})*xnew{j};
        %sd = sd + dot(a1,a2);
        %sn = sn + norm(a2)^2;
        sd = dot(a1,a2);
        sn = norm(a2)^2;
    if sn~=0
        c2 = sd/sn;
        x{j}=x{j}+xnew{j}*c2;
        if minc
            x{j}=max(x{j},grays(1));
        end
    end
end
end
    
end     

