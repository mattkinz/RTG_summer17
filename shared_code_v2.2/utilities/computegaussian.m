function T = computegaussian(rx,rz,sigmax,sigmaz,slices,chunksize)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


if slices<rz || chunksize<rz
    rz=1;
end

if rz == 0 
    rz=1;
end
T = zeros(rx,rx,rz);
middle = [rx/2+.5, rz/2+.5];
for i = 1:rx
    for j  = 1:rx
        for k = 1:rz
            val = ((middle(1)-i)^2+(middle(1)-j)^2)/sigmax+(middle(2)-k)^2/sigmaz;
            T(i,j,k) = exp(-val);
        end
    end
end

T = T/sum(sum(sum(T)));

end