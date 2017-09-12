function x = sirt_sub_chunk(W,Wn,bb,iter,x,set,c,minc,grays)


% Written by: Toby Sanders @ Pacific Northwest National Laboratory
% Computational & Applied Mathematics Department, Univ. of South Carolina
% 7/11/2014



slices = max(size(set));
xnew = cell(slices,1);

for i = 1:slices
    subdim=max(size(set{i}));
    xnew{i}=zeros(subdim,1);
end

for i = 1:iter
    sn=0;sd=0;
    for j = 1:slices       
        a1=bb(:,j)-W(:,set{j})*x{j};
        %xnew{j}=(a1'*Wn(:,set{j}))'.*c{j};
        xnew{j} = (W(:,set{j})'*(Wn.*a1)).*c(set{j});
        a2=W(:,set{j})*xnew{j};
        sd = sd + dot(a1,a2);
        sn = sn + norm(a2)^2;
    end
    if sn~=0
        c2 = sd/sn;
    for j = 1:slices
        x{j}=x{j}+xnew{j}*c2;
        if minc
            x{j}=max(x{j},grays(1));
        end
    end    
    end
end
    
end     

