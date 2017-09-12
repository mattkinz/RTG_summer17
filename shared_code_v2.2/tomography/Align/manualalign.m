function [stack,shifts] = manualalign(stack)

%MANUALALIGN%

%DESCRIPTION:  An interactive algorithm for manual alignment.

%NOTATION: [stacknew,s] = manualalign(stack);

%INPUT: 
    %stack - tilt series
    
%OUTPUTS:  
    %stacknew - a new stack, with shifts based on user inputs
    %s - the shifts that were applied
    
%Notes:
%This function is for interactive alignment.
%The stack of images will be displayed in pairs, and will play as a movie.
%After a pair is shown, the movie pauses and the user may input 
%a command for shifting the image shown.  If an input is not given, the 
%movie will simply continue to the next pair.
%The inputs for shifting are "u" for up, "d" for down, "l" for left, 
%and "r" for right

%Once a command is entered for a shift, more shift inputs will be taken
%and the movie will not continue until the user either inputs "apply"
%to apply the shifts or "cancel", to not apply the shifts.
%If the user inputs "apply", then the user will be asked if the shifts 
%should be applied to all of the proceeding (or preceeding, depending) images
%Once the user is satisfied of the alignment, input "finish"
%to close the algorithm and return the aligned stack.


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


s =0;
direction='f';
move='';
k = size(stack,3);
temp = imresize(stack,[512,512]);
shifts = zeros(k,2);
fprintf('Press enter to scroll through the image pairs\n');
fprintf('To switch the direction of the scrolling press:\n');
fprintf('"b" to scroll back\n');
fprintf('"f" to scroll forward\n\n');
fprintf('To begin shifting an image press:\n')
fprintf('"u" for up\n')
fprintf('"d" for down\n')
fprintf('"l" for left\n')
fprintf('"r" for right\n\n')
fprintf('Once satisfied with the alignment, input "finish"\n');
fprintf('to stop the operation and return the new stack\n');

while ~strcmp(move,'finish')      
    if strcmp(direction,'f')
        s=mod(s+1,k);
        if s==0
            s=k;
            direction='b';
        end
    elseif strcmp(direction,'b')
        s=mod(s-1,k);
        if s==1
            direction='f';
        end
    end
    if strcmp(direction,'f')
        figure(1);
        imshowpair(temp(:,:,s),temp(:,:,s+1));
        title(sprintf('images %i and %i',s,s+1));
    else
        figure(1);
        imshowpair(temp(:,:,s),temp(:,:,s-1));
        title(sprintf('images %i and %i',s,s-1));
    end
    
    move = input('','s');
    if sum(strcmp(move,{'u','d','l','r'}))
        if strcmp(direction,'f')
            [shift,apply]=pairmove(temp(:,:,s),temp(:,:,s+1),move,s,s+1);
            if strcmp(apply,'apply')
                applyall=input('Apply the shifts to all of the proceeding images? y/n:','s');
                if strcmp(applyall,'y')                
                    for i = s+1:k
                        shifts(i,:)=shifts(i,:)+shift;
                        temp(:,:,i)=circshift(temp(:,:,i),shift);
                    end
                    fprintf('Shifts applied to all of the proceeding images\n');
                else
                    shifts(s+1,:)=shifts(s+1,:)+shift;
                    temp(:,:,s+1)=circshift(temp(:,:,s+1),shift);
                end
            end
        else
            [shift,apply]=pairmove(temp(:,:,s),temp(:,:,s-1),move,s,s-1);
            if strcmp(apply,'apply')
                applyall=input('Apply the shifts to all of the preceeding images? y/n:','s');
                if strcmp(applyall,'y') 
                    for i = s-1:1
                        shifts(i,:)=shifts(i,:)+shift;
                        temp(:,:,i)=circshift(temp(:,:,i),shift);
                    end
                    fprintf('Shifts applied to all of the preceeding images\n');
                else
                    shifts(s-1,:)=shifts(s-1,:)+shift;
                    temp(:,:,s-1)=circshift(temp(:,:,s-1),shift);
                end
            end
        end
    elseif sum(strcmp(move,{'f','b'}))
        direction=move;
    end
    
            
                     
end

shifts(:,1)=shifts(:,1)*size(stack,1)/512;
shifts(:,2)=shifts(:,2)*size(stack,2)/512;
shifts = round(shifts);
for i =1:k
    stack(:,:,i)=circshift(stack(:,:,i),shifts(i,:));
end


end
