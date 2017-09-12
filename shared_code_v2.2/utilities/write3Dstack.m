function [] = write3Dstack(rec,fname,dim)



%WRITE_STACK%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DESCRIPTION:
    %This function is for writing a 3D image file for a 3D volume.
    
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 04/18/2016

if ~exist('fname','var')
    fname = 'imgfile.tif';
end

if ~exist('dim','var')
    dim = 3;
end

k = size(rec,dim);

for i = 1:k
    if dim==3
        imwrite(rec(:,:,i),fname,'WriteMode','append');
    elseif dim==2
        imwrite(squeeze(rec(:,i,:)),fname,'WriteMode','append');
    else
        imwrite(squeeze(rec(i,:,:)),fname,'WriteMode','append');
    end
end


end