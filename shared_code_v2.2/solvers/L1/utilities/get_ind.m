function ind = get_ind(k,p,q,r)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


U = zeros(p,q,r);
if k<=q
    U(:,end-k+1:end,:)=1;
    indx = find(U);
    U(:)=0;
else
    indx=[];
end
if k<=p
    U(end-k+1:end,:,:)=1;
    indy = find(U);
    U(:)=0;
else
    indy=[];
end
if k<=r
    U(:,:,end-k+1:end)=1;
    indz = find(U);
else
    indz=[];
end

ind = [indx;indy+p*q*r;indz+2*p*q*r];