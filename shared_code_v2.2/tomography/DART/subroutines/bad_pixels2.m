function [set,fixed] = bad_pixels2(W,bb,U,tol)

%function for determining the voxels involved  with the measurements with
%large error


[m,n,k]=size(U);
set = cell(k,1);
fixed = cell(k,1);
error = (W*reshape(U,m*n,k)-bb);
Tr = computegaussian(5,5,1,3,k,k);

[numray,k]=size(bb);
numget = round(m*n*tol);
v = W'*error;
%figure(10);imagesc(reshape(v,m,n));figure(10);
v = abs(v);
v = reshape(v,m,n,k);
v = imfilter(v,Tr);
v = reshape(v,m*n,k);
U = reshape(U,m,n,k);
mask = zeros(m,n,k,'single');
f = find([diff(U,1,2),zeros(m,1,k)]);
mask(f)=1;
mask(f+m)=1;
f = find([diff(U,1,1);zeros(1,n,k)]);
mask(f)=1;
mask(f+1)=1;
f = find(cat(3,diff(U,1,3),zeros(m,n)));
mask(f)=1;
mask(f+m*n)=1;
mask = imfilter(mask+.1,Tr);
v = v.*reshape(mask,m*n,k);

for i =1:k
    [s1,s2]=sort(v(:,i));
    s2 = wrev(s2);
    set{i} = s2(1:numget);
    fixed{i} = s2(numget+1:end);
end

end
