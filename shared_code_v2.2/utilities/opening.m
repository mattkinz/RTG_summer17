function T = opening(volume,object)

%This function is for erosion and dialation procedure.

%volume should be binary
%object should be binary as well, but will work either way

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


[m,n,p]=size(volume);
%volume = logical(volume);
T = zeros(m,n,p);
tic;
F = imfilter(volume,object);
toc;
tol = sum(sum(sum(object)));
[a,b,c]=size(object);
nz = find(object);
mm = min(object(nz));
tol2 = tol-mm;

F = max(F-tol2,0);

T = threshold(F,.001,[0 1]);
%{
for i = 1:max(size(nz))
    k0 = floor(nz(i)/a/b)+1;
    j0 = floor((nz(i)-(k0-1)*a*b)/a)+1;
    i0 = nz(i)-(j0-1)*a-(k0-1)*a*b;
    nz(i)=i0+(j0-1)*m+(k0-1)*m*n;
end
ra= floor(a/2);
rb = floor(b/2);
rc = floor(c/2);

whos
for i = 1:m-a+1
    for j = 1:n-b+1
        for k = 1:p-c+1
            if F(i,j,k)
                T((i-ra)+(j-1-rb)*m+(k-1-rc)*m*n+nz)=1;
            end
        end
    end
end

end
%}