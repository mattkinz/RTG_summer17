function [set,fixed] = bad_pixels1(W,bb,U,tol,angles)

%function for determining the voxels involved  with the measurements with
%large error

angles = angles*pi/180;


[m,n,k]=size(U);
set = cell(k,1);
fixed = cell(k,1);
error = (W*reshape(U,m*n,k)-bb)./(bb+mean(bb));
Tr = computegaussian(3,3,3,3,k,k);

[numray,k]=size(bb);
numray = numray/max(size(angles));
for i =1:max(size(angles))
    error((i-1)*numray+1:i*numray)=error((i-1)*numray+1:i*numray)*abs(sin(angles(i)));
end

numget = round(m*n*tol);
v = W'*abs(error);
figure(10);imagesc(reshape(v,m,n));figure(10);
v = abs(v);
v = reshape(v,m,n,k);
v = imfilter(v,Tr);
v = reshape(v,m*n,k);
for i =1:k
    [s1,s2]=sort(v(:,i));
    s2 = wrev(s2);
    set{i} = s2(1:numget);
    fixed{i} = s2(numget+1:end);
end

end
