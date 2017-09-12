function x = sirt_sub_chunk_reg(W,M,bb,iter,x,set,fixed,minc,grays,lambda1,lambda2)


%Written by: Toby Sanders @
%Computational & Applied Mathematics Department, Univ. of South Carolina
%6/9/15


slices = max(size(set));
xnew = cell(slices,1);
xf = cell(slices,1);
subdim=zeros(slices,1);
for i = 1:slices
    subdim(i)=max(size(set{i}));
    xnew{i}=zeros(size(W,2),1);
    xf{i}=x(fixed{i},i);
end


re1 = @(x)reshape(x,size(x,1)*size(x,2)*size(x,3),1);
re2 = @(x)reshape(x,sqrt(size(x,1)),sqrt(size(x,1)));
%lambda1=5;lambda2=10;


for i = 1:iter
    sn=0;sd=0;
    for j = 1:slices       

        z1 = bb(:,j)-W(:,set{j})*x(set{j},j);
        z2 = -lambda1*M*x(:,j);
        z3 = lambda2*(xf{j}-x(fixed{j},j));
        
        
        xnew{j}(set{j})=W(:,set{j})'*z1;
        xnew{j}(fixed{j})=z3*lambda2;
        xnew{j}=xnew{j}+M'*z2*lambda1;
            
        y1=W(:,set{j})*xnew{j}(set{j});
        y3=xnew{j}(fixed{j})*lambda2;
        y2=lambda1*M*xnew{j};

        sd = sd + dot([z1;z2;z3],[y1;y2;y3]);
        sn = sn + norm([y1;y2;y3])^2;

    end
    if sn~=0
        c2 = sd/sn;
        for j = 1:slices
            x(:,j)=x(:,j)+xnew{j}*c2;
            if minc
                x(:,j)=max(x(:,j),grays(1));
            end
            xnew{j}(:)=0;
        end    
    end
end
    
end     

