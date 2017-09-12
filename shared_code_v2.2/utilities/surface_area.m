function [surfacearea,mask] = surface_area(V)


%SURFACE_AREA%%%%%%%%%%%%%%%%%%%%%%


%DESCRIPTION:
       %This function computes the surface area of binary 3D volume.
       
%NOTATION:
    %[sa,mask] = surface_area(V);
    
%INPUTS:
    %V - 3 dimensional volume matrix with ONLY 2 gray values.
    
%OUTPUTS:
    %sa - the total number of boundary faces in V
    %mask - a mask of the boundary
    
%NOTES:
    %The output sa is simply sum(sum(mask)).  mask should actually have
    %several gray values, as some voxels will have multiple boundary faces.
    %Make sure that V has a smooth surface, or there will be far too many
    %boundary faces in mask.  To make sure that this is so, look at the
    %mask and make sure the mask has think boundaries.

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


[m,n,k]=size(V);

mask = zeros(m*n*k,1,'uint8');

f = find(reshape([diff(V,1,2),zeros(m,1,k)],m*n*k,1));
mask(f)=mask(f)+1;
f = find(reshape([diff(V,1,1);zeros(1,n,k)],m*n*k,1));
mask(f)=mask(f)+1;
f = find(reshape(cat(3,diff(V,1,3),zeros(m,n)),m*n*k,1));
mask(f)=mask(f)+1;

mask = reshape(mask,m,n,k);

surfacearea = sum(sum(sum(mask)));

end

