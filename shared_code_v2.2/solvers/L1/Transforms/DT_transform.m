function [D,Dt] = DT_transform(k,p,q,r,init,g)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

scl = 1;
V = get_DT_weights(init,g);
D = @(U)D_Forward(U,k,p,q,r,V,scl);
Dt = @(Uc)D_Adjoint(Uc,k,p,q,r,V,scl);


% high order finite differences
function [dU] = D_Forward(U,k,p,q,r,V,scl)


U = reshape(U,p,q,r);
if k~=0
    if k<=size(U,2)
        X = diff(U,k,2);
        x = zeros(size(U,1),k,size(U,3));
        for i = 1:k
            x(:,i,:) = diff([U(:,end-k+i:end,:),U(:,1:i,:)],k,2);
        end
        X = [X,x];
    else
        X = zeros(size(U));
    end

    if k<=size(U,1)
        Y = diff(U,k,1);
        y = zeros(k,size(U,2),size(U,3));
        for i = 1:k
            y(i,:,:) = diff([U(end-k+i:end,:,:);U(1:i,:,:)],k,1);
        end
        Y = [Y;y];
    else
        Y = zeros(size(U));
    end


    if k<=size(U,3)
        Z = diff(U,k,3);
        z = zeros(size(U,1),size(U,2),k);
        for i = 1:k
            temp = cat(3,U(:,:,end-k+i:end),U(:,:,1:i));
            z(:,:,i) = diff(temp,k,3);
        end
        Z = cat(3,Z,z);
    else
        Z = zeros(size(U));
    end
else
    %standard l1 minimization of order 0
    X=U;Y=U;Z=U;
end
ng = size(V,2);
dU = zeros(p*q*r,3+ng);
dU(:,1:3) = [X(:),Y(:),Z(:)];
U = U(:);


V = scl*V;
for i = 1:ng
    dU(:,3+i) = U(:).*V(:,i);
end


%transpose FD
function dtxy = D_Adjoint(dU,k,p,q,r,V,scl)
X = reshape(dU(:,1),p,q,r);
Y = reshape(dU(:,2),p,q,r);
Z = reshape(dU(:,3),p,q,r);

ng = size(V,2);

U = zeros(p*q*r,ng);
V = scl*V;
for i = 1:ng
    U(:,i) = dU(:,3+i).*conj(V(:,i));
end

if k~=0
    if k<=size(X,2)
        Dtx = (-1)^k*diff(X,k,2);
        dtx = zeros(size(X,1),k,size(X,3));
        for i = 1:k
            dtx(:,i,:)= (-1)^k*diff([X(:,end-k+i:end,:),X(:,1:i,:)],k,2);
        end
        Dtx = [dtx,Dtx];
    else
        Dtx = 0;
    end



    if k<=size(Y,1)
        Dty = (-1)^k*diff(Y,k,1);
        dty = zeros(k,size(Y,2),size(Y,3));
        for i = 1:k
            dty(i,:,:) = (-1)^k*diff([Y(end-k+i:end,:,:);Y(1:i,:,:)],k,1);
        end

        Dty = [dty;Dty];
    else
        Dty = 0;
    end


    if k<=size(Z,3)
        Dtz = (-1)^k*diff(Z,k,3);
        dtz = zeros(size(Z,1),size(Z,2),k);
        for i =1:k
            dtz(:,:,i) = (-1)^k*diff(cat(3,Z(:,:,end-k+i:end),Z(:,:,1:i)),k,3);
        end

        Dtz = cat(3,dtz,Dtz);
    else
        Dtz = 0;
    end

    
    dtxy = Dty + Dtx + Dtz;
    dtxy = dtxy(:) + sum(U,2);
else
    %standard l1 minimization
    dtxy = X(:)+Y(:)+Z(:) + sum(U,2);
end

