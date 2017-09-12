function Gamma = Gamma_matrix(angles,slices)

k = max(size(angles));
Gamma = sparse(k*slices,k*slices);

a = angles*pi/180;
C1 = sum(cos(a).^2);
S1 = sum(sin(a).^2);
T = -sum(sin(a).*cos(a));
A = [C1,T;T,S1];
Phi = zeros(k,2);
Phi(:,1)=cos(a);
Phi(:,2)=-sin(a);

G = Phi*inv(A)*Phi'-eye(k);

n = slices;
y = zeros(k^2*n,3);

for i =1:k
    for j =1:n
        y((i-1)*k+(j-1)*k^2+1:i*k+(j-1)*k^2,1)=i+(j-1)*k;
        y((i-1)*k+(j-1)*k^2+1:i*k+(j-1)*k^2,2)=(j-1)*k+1:j*k;
        y((i-1)*k+(j-1)*k^2+1:i*k+(j-1)*k^2,3)=G(i,:);
    end
end

Gamma = sparse(y(:,1),y(:,2),y(:,3));

end