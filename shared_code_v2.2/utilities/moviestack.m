function [] = moviestack(stacks,cycles,direction,map)


%MOVIESTACK and MOVIESTACKG

%DESCRIPTION:
    %these functions are for 2D visualizations of 3D volumes.  They will
    %scroll through orthogonal slices of an input 3D volume.
    
%NOTATION:
    %moviestack(stack,cycles,direction,map)
    
%INPUTS:
    %stack - 3D volume to be visualized
    %cycles - number of cycles through the volume.  Default is 10.
    %direction - specifies the direction which to scroll through the
        %volume.  Set to either 1, 2, or 3.  1 will scroll through the z
        %direction, 2 will scroll through the x direction and 3 will scroll
        %through the y direction.  Default is 1.
    %map - a string input specifying the colormap.  Some of the best are
    %'gray', 'jet', 'hot', and 'bone'.  Default is 'gray'.
    
    
%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

%last update: 02/26/2016


if ~exist('direction','var')
    direction=1;
end
if ~exist('cycles','var')
    cycles=10;
end
if ~exist('map','var')
    map = 'gray';
elseif ~ischar(map)
    map = 'gray';
elseif ~sum(strcmp(map,{'jet','gray','hot','bone','cool','hsv','spring',...
        'summer','autumn','winter','copper','pink','lines','parula'}))
    map = 'gray';
end


if isa(stacks,'cell')
    if direction ==2
       for i = 1:cycles
        if mod(i,2)==1
            for j = 1:size(stacks{1},2)
                for k = 1:max(size(stacks))
                    figure(k);imagesc(squeeze(stacks{k}(:,j,:)));
                    colormap(map);colorbar;
                    title(j);
                end
                pause;
            end
        else
            for j = size(stacks{1},2):-1:1
                for k = 1:max(size(stacks))
                    figure(k);imagesc(squeeze(stacks{k}(:,j,:)));
                    colormap(map);colorbar;
                    title(j);
                end
                pause;
            end
        end
      end
    elseif direction ==3
        for i = 1:cycles
        if mod(i,2)==1
            for j = 1:size(stacks{1},1)
                for k = 1:max(size(stacks))
                    figure(k);imagesc(squeeze(stacks{k}(j,:,:)));
                    colormap(map);colorbar;
                    title(j);
                end
                pause;
            end
        else
            for j = size(stacks{1},1):-1:1
                for k = 1:max(size(stacks))
                    figure(k);imagesc(squeeze(stacks{k}(j,:,:)));
                    colormap(map);colorbar;
                    title(j);
                end
                pause;
            end
        end
        end
    else
        for i = 1:cycles
            if mod(i,2)==1
                for j = 1:size(stacks{1},3)
                    for k = 1:max(size(stacks))
                        figure(k);imagesc(stacks{k}(:,:,j));
                        colormap(map);colorbar;
                        title(j);
                    end
                    pause;
                end
            else
                for j = size(stacks{1},3):-1:1
                    for k = 1:max(size(stacks))
                        figure(k);imagesc(stacks{k}(:,:,j));
                        colormap(map);colorbar;
                        title(j);
                    end
                    pause;
                end
            end
        end
    end
else
    if direction ==2
        for i = 1:cycles
            if mod(i,2)==1
                for j = 1:size(stacks,2)
                    figure(1);imagesc(squeeze(stacks(:,j,:)));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            else
                for j = size(stacks,2):-1:1
                    figure(1);imagesc(squeeze(stacks(:,j,:)));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            end
        end
    elseif direction ==3
        for i = 1:cycles
            if mod(i,2)==1
                for j = 1:size(stacks,1)
                    figure(1);imagesc(squeeze(stacks(j,:,:)));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            else
                for j = size(stacks,1):-1:1
                    figure(1);imagesc(squeeze(stacks(j,:,:)));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            end
         end
    else
        for i = 1:cycles
            if mod(i,2)==1
                for j = 1:size(stacks,3)
                    figure(1);imagesc(stacks(:,:,j));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            else
                for j = size(stacks,3):-1:1
                    figure(1);imagesc(stacks(:,:,j));
                    colormap(map);colorbar;
                    title(j);pause;
                end
            end
        end
    end
end