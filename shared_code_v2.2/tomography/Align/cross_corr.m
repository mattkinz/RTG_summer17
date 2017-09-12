function [stack,shifts] = cross_corr(stack,startslice)


%DESCRIPTION
    %This function performs cross correlation of an input stack of images.
    %
%NOTATION
    %[stacknew,shifts]=cross_corr(stack,startslice);
    %
%INPUTS
    %stack - 3D matrix assumed to be a stack of images to be cross 
        %correlated
    %startslice - the reference image that all slices are consecutively
    %aligned to.  This image will not move
    %
%OUTPUTS
    %stack- 3D stack of images after cross correlation
    %shifts - shifts that were applied to each image
    %
%NOTE: This function automatically uses a window to improve alignment.  To
    %remove the windowing, set pwr=0.
   %
%
%   
%
% Written by: Toby Sanders @ Pacific Northwest National Laboratory
% Computational & Applied Mathematics Department, Univ. of South Carolina
% 7/11/2014



[m,n,k]=size(stack);
if ~exist('startslice','var')
    startslice=ceil(k/2);
end
shifts = zeros(k,2);
w = zeros(m,n);
a = (m+1)/2;
b = (n+1)/2;
pwr = 0;
if pwr~=0
    fprintf('Note that windowing is being used!!! \n');
end
for y =1:m
    for x = 1:n
        w(y,x)=(1-(y-a)^2/2/a^2-(x-b)^2/2/b^2)^(pwr);
    end
end

for i = startslice+1:k
    t1 = fft2(stack(:,:,i-1).*w);
    t2 = fft2(stack(:,:,i).*w);
    cc = ifft2(t1.*conj(t2));
    [~,b]=max(cc);
    [~,c]=max(max(cc));
    s1 = b(c)-1;
    s2 = c-1;
    shifts(i,1)=mod(s1+shifts(i-1,1),m);
    shifts(i,2)=mod(s2+shifts(i-1,2),n);
end

for i = startslice-1:-1:1
    t1 = fft2(stack(:,:,i).*w);
    t2 = fft2(stack(:,:,i+1).*w);
    cc = ifft2(t1.*conj(t2));
    [~,b]=max(cc);
    [~,c]=max(max(cc));
    s1 = -b(c)+1;
    s2 = -c+1;
    shifts(i,1)=mod(s1+shifts(i+1,1),m);
    shifts(i,2)=mod(s2+shifts(i+1,2),n);
end

% mod the shifts
bigg1 = find(shifts(:,1)>m/2);
bigg2 = find(shifts(:,2)>n/2);
shifts(bigg1,1) = shifts(bigg1,1)-m;
shifts(bigg2,2) = shifts(bigg2,2)-n;

% apply shifts
for i = 1:k
    stack(:,:,i)=circshift(stack(:,:,i),shifts(i,:));
end



end