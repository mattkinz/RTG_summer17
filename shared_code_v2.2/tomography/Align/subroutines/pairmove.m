function [shift,move] = pairmove(image1,image2,move,s1,s2)

%subroutine for manualalign


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


fprintf('Beginning shifting mode for image %i overlayed with image %i\n',s2,s1);
fprintf('after shifting is complete, input "apply" to apply the shifts\n')
fprintf('or "cancel" to cancel the shifts\n');
shift = zeros(1,2);
while ~sum(strcmp(move,{'apply','cancel'}))
    
    move = strcmp(move,{'u','d','l','r','apply','cancel'}).*...
        [-1 1 -1 1 1 1];
    if sum(move)==0
        fprintf('Invalid input\n');
    end
    move = [sum(move(1:2)),sum(move(3:4))];
    image2 = circshift(image2,move);
    imshowpair(image1,image2);
    title(sprintf('images %i and %i',s1,s2));
    shift = shift+move;
    move = input('','s');
end

