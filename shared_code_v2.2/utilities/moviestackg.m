function [] = moviestackg(stacks,cycles,direction)


%MOVIESTACK and MOVIESTACKG

%DESCRIPTION:
    %these functions are for 2D visualizations of 3D volumes.  They will
    %scroll through orthogonal slices of an input 3D volume.
    
%NOTATION:
    %moviestackg(stack,cycles,direction)
    
%INPUTS:
    %stack - 3D volume to be visualized
    %cycles - number of cycles through the volume.  Default is 10.
    %direction - specifies the direction which to scroll through the
        %volume.  Set to either 1, 2, or 3.  1 will scroll through the z
        %direction, 2 will scroll through the x direction and 3 will scroll
        %through the y direction.  Default is 1.

    
    
%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



if ~exist('direction','var')
    direction=1;
end
if ~exist('cycles','var')
    cycles=10;
    fprintf('Number of cycles not specified, cycling for 10\n');
end

if isa(stacks,'cell')
    if direction ==2
        for i = 1:max(size(stacks))
            [m,n,k]=size(stacks{i});
            temp = zeros(m,k,n);
            for j =1:n
                temp(:,:,j)=stacks{i}(:,j,:);
            end
            stacks{i}=temp;
        end
        clear temp;
    elseif direction ==3
        for i = 1:max(size(stacks));
            [m,n,k]=size(stacks{i});
            temp = zeros(n,k,m);
            for j = 1:m
                temp(:,:,j)=stacks{i}(j,:,:);
            end
            stacks{i}=temp;
        end
        clear temp;
    end
    for i = 1:cycles
        if mod(i,2)==1
            for j = 1:size(stacks{1},3)
                for k = 1:max(size(stacks))
                    figure(k);imshow(stacks{k}(:,:,j));
                    title(j);
                end
                pause;
            end
        else
            for j = size(stacks{1},3):-1:1
                for k = 1:max(size(stacks))
                    figure(k);imshow(stacks{k}(:,:,j));
                    title(j);
                end
                pause;
            end
        end
    end
else
    if direction ==2
        [m,n,k]=size(stacks);
        temp = zeros(m,k,n);
        for i =1:n
            temp(:,:,i)=stacks(:,i,:);
        end
        stacks=temp;
        clear temp;
    elseif direction ==3
        [m,n,k]=size(stacks);
        temp = zeros(n,k,m);
        for i = 1:m
            temp(:,:,i)=stacks(i,:,:);
        end
        stacks=temp;
        clear temp;
    end
    for i = 1:cycles
        if mod(i,2)==1
            for j = 1:size(stacks,3)
                figure(1);imshow(stacks(:,:,j));
                title(j);pause;
            end
        else
            for j = size(stacks,3):-1:1
                figure(1);imshow(stacks(:,:,j));
                title(j);pause;
            end
        end
    end
end
end