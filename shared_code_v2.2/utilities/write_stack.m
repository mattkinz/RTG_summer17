function [] = write_stack(rec,type)



%WRITE_STACK%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DESCRIPTION:
    %This function is for writing images of a 3D volume.  The images will
    %be written in the current open folder in MATLAB.
    
%NOTATION:
    %write_stack(stack,type)
    %write_stack(rec,type)
    
%INPUTS:
    %rec, or stack - the 3D volume or stack of images that will be written
        %out of MATLAB as image files
    %type - a character specifying the image file type that will be
        %written.  Examples include 'bmp','jpg','jpeg','pbm','png', 'tif'.
        %Default is 'jpg'.
%OUPUTS:
    %No outputs, you will just have a bunch of images now in your current
    %folder.
    
    
%NOTES:
    %The MATLAB function "imwrite" is used for this function.  This
    %function writes images in the same way that "imshow" shows images.  So
    %make sure that the image intensities are scaled accordingly.
    

    
%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

if ~exist('type','var')
    type='png';
elseif ~ischar(type)
    type='png';
elseif ~sum(strcmp(type,{'bmp','jpg','jpeg','pbm','png',...
            'tif'}))
        type='jpg';
end

fprintf(['\nWriting ',type,' images\n\n']);

k=size(rec,3);

for i =1:k
    if i<10
        imwrite(rec(:,:,i),sprintf(['000',num2str(i),'.',type]));
    elseif i<100
        imwrite(rec(:,:,i),sprintf(['00',num2str(i),'.',type]));
    elseif i<1000
        imwrite(rec(:,:,i),sprintf(['0',num2str(i),'.',type]));
    else
        imwrite(rec(:,:,i),sprintf([num2str(i),'.',type]));
    end
end

end