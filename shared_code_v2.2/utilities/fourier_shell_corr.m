function fsc = fourier_shell_corr(v1,v2,epsilon,voxel_size)


%DESCRIPTION: 
    %This function computes the N-dimensional Fourier shell correlation of
    %two N-dimensional volumes.
    
%NOTATION:
    %fsc = fourier_shell_corr(v1,v2,epsilon,voxel_size);

%INPUTS: 
    %v1,v2 - input volumes for computing the FSC.
    %epsilon - thickness (in number of voxels or frequency levels) of the 
        %shells in the fourier domain.  This variable is optional.  
        %Default value is 3.
    %voxel_size - length of the edges of the voxels in nanometers.  This is
        %just used for the plot, and is optional.  Default value is 1.
    
%OUTPUT: 
    %fsc - the fourier shell correlation of the two volumes, along with a
        %plot of the fsc.
        
        

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014
        
        
if ~exist('voxel_size','var')
    voxel_size=1;
end
if ~exist('epsilon','var')
    epsilon=3;
end

[d1,d2,d3]=size(v1);
[d4,d5,d6]=size(v2);
if d1~=d4 || d2~=d5 || d3~=d6
    error('volume dimensions do not agree');
end
if d3==1
    n=2;
    if d2==1
        n=1;
    end
else
    n=3;
end
fprintf('Computing Fourier transforms\n\n');
f1 = fftshift(fftn(v1));
f2 = fftshift(fftn(v2));
f3 = real(f1.*conj(f2));


levels = floor(max([d1,d2,d3])/epsilon/2);
n1 = zeros(levels ,1);
n2 = zeros(levels ,1);
fsc = zeros(levels ,1);
fprintf('Computing the FSC\n\n');
for j = 0:d1-1
    d1t = ((d1+1)/2-j-1/2)^2;
    for k = 0:d2-1
        d2t = d1t+((d2+1)/2-k-1/2)^2;
        for l = 0:d3-1
            dist = sqrt(d2t+((d3-1)/2-l-1/2)^2);
            level = floor(dist/epsilon)+1;
            if level < levels +1
                fsc(level)=f1(j+1,k+1,l+1).*conj(f2(j+1,k+1,l+1))+fsc(level);
                %fsc(level)=f3(j+1,k+1,l+1)+fsc(level);
                n1(level) = n1(level)+abs(f1(j+1,k+1,l+1))^2;
                n2(level) = n2(level)+abs(f2(j+1,k+1,l+1))^2;        
            end
        end
    end
end



fsc=fsc./sqrt(n1);
fsc=fsc./sqrt(n2);

v = 0:1/voxel_size/(max(size(fsc))-1):1/voxel_size;
%whos
%plot(v,real(fsc))
%xlabel('1/nm')
%ylabel('FSC')

end
    
        

