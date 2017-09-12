function [D,Dt] = FD3D_dips(k,p,q,r)


% Written by Toby Sanders @ASU
% Department of Statistical & Mathematical Sciences
% 05/17/2016


% finite difference operators for polynomial annihilation
% k is the order of the PA transform
D = @(U)D_Forward(U,k,p,q,r);
Dt = @(Uc)D_Adjoint(Uc,k,p,q,r);


% high order finite differences
function [dU] = D_Forward(U,k,p,q,r)

if k~=0
    U = reshape(U,p,q,r);
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
    dU = [X(:),Y(:),Z(:)];
else
    %standard l1 minimization of order 0
    dU = U(:);
end



%transpose FD
function dtxy = D_Adjoint(dU,k,p,q,r)

if k~=0
    X = reshape(dU(:,1),p,q,r);
    Y = reshape(dU(:,2),p,q,r);
    Z = reshape(dU(:,3),p,q,r);
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
    dtxy = dtxy(:);%*2^(1-k);
else
    %standard l1 minimization
    dtxy = dU(:);
end
    
    
    
    