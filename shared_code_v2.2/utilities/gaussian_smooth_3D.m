function V = gaussian_smooth_3D(V,rx,rz,sigmax,sigmaz)


%GAUSSIAN_SMOOTH_3D%%%%%%%%%%%%%%%%%%%%%%

%DESCRIPTION:
    %This function is for 3D smoothing of a volume.  It also works for 2D
    %images.
    
%NOTATION:
    %Vnew = guassian_smooth_3D(V,rx,rz,sigmax,sigmaz);
    
%INPUTS:
    %V - 3D volume to be smoothed
    %rx - radius of the gaussian function for averaging in the x and y
    %coordinates.
    %rz - radius of the averaging function for averaging in the z
    %coordinate.
    %sigmax - sigma value for x and y of the gaussian function.
    %sigmaz - sigma value for z of the gaussian function.
    
%NOTES:
    %rx and rz set to 3 or 5 are usually sufficient.  sigmax and sigmaz set
    %to 1 or 2 is usually sufficient.  It all depends on the noise in the
    %volume.  Higher values for sigmax and sigmaz will create larger
    %smoothing.  For a reconstructed tomography volume, there may be
    %inconsistencies in the z-direction, therefore in this case it can be
    %better to set sigmaz higher than sigmax.  There are many ways to
    %filter images using MATLAB.  This function is just a direct convenient
    %filter that should work well for most tomography reconstructions.
    
%OUTPUTS:
    %Vnew - the new smoothed function.

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014
    
    
    Tr = computegaussian(rx,rz,sigmax,sigmaz,size(V,3),size(V,3));

    V = imfilter(V,Tr);

end