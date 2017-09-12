function [stack,shifts,usables] = COM_align(stack,angles,ratio,s)


%  Center of mass alignment for electron tomography data
%
% function [stack,shifts,usables] = COM_align(stack,angles,ratio,s)
%
% Inputs: stack is 3-D matrix containing tomography data 
% angles are the angles in degrees
% ratio is the ratio of slices to be used for each projection
% s is the number of projections used as a sequence to determine rigid
% alignment.  Setting s to be the number of angles is equivalent to the
% alignment method published here: 
% http://ascimaging.springeropen.com/articles/10.1186/s40679-015-0005-7
%
%
%
% The ratio variable:
% should be input as a number between 0 and 1, where the user rates the
% consitency of the stacks through the slices.  From our experience, very 
% small values of s are better.  If the user does not input a
% rating, the default value of .1 is used.  Setting ratio
% to perfect 1 would mean that through the slices, there is significant amount 
% of mass in all of the slices, and very little mass moves in or out as the
% tilt angle changes.  A lower ratio would mean otherwise.
% This rating is then used to determine how many slices to use for the
% center of mass alignment.  For example, if you set the rating to .2, then
% only the best 20 percent of the slices will be used in the alignment step.
% That way, if the the stack has some slices with very little mass, yet 
% mass moves into these slices at high angles, these slices may be ignored 
% for the aligment.
%
%
% Written by Toby Sanders @ASU & @PNNL & @USC
% School of Mathematical & Statistical Sciences
% 06/08/2016



% Assign default values to unspecified parameters
if (nargin < 3 || isempty(ratio)), ratio = 1/10; end
if (nargin < 4 || isempty(s)), s = numel(angles); end;



%check the dimensions and angles, and convert to radians
[m,n,k0] = size(stack);
[k,ind] = max(size(angles));
if ind ==2,angles=angles';end
if k~=k0
    error('size of stack doesnt match input angles');
end    
angles = angles*pi/180;




%store the stack into a cell to align with the sequences
stacks = cell(k-s+1,1);
for i = 1:k-s+1
    stacks{i} = stack(:,:,i:i+s-1);
end


%initialize vectors that will hold COM computations
%variable t holds the COM location
%variable ss holds the total mass
t=zeros(s,k-s+1,n);
ss = zeros(s,k-s+1,n);


%variable w is for the COM calculation
w = (1:m)';

%Compute all centers of mass for each sequence
for i = 1:n
    for j = 1:k-s+1
        for l = 1:s
            ss(l,j,i) = sum(stacks{j}(:,i,l));
            t(l,j,i) = sum(stacks{j}(:,i,l).*w)/ss(l,j,i);
        end
    end
end
%shift so that the center is 0
t = t-(m+1)/2;



%compute mean absolute difference divided by the median mass
ss2 = median(ss);
for l = 1:s
    ss(l,:,:) = abs((ss(l,:,:)-ss2)./ss2);
end
ss2 = mean(ss);


%determine the best slices for each sequence
%store the COM location of these into t_select
%The "best" slices those in which ss2 is minimal, 
%i.e. mean absolute difference divided by the mass is minimal
num = round(ratio*n);
usables = zeros(num,k-s+1);
t_select = zeros(s*num,(k-s+1));
disp_mat = zeros(n,k-s+1);
for i = 1:k-s+1
    [~,s3]=sort(ss2(1,i,:));
    usables(:,i) = reshape(s3(1:num),num,1);
    t_select(:,i)=reshape(t(:,i,usables(:,i)),s*num,1);
    disp_mat(usables(:,i),i)=1;
end



%display the mean absolute difference matrix, and the slices chosen for
%alignment of each sequence

figure(27);

subplot(2,1,1);
imagesc(squeeze(log10(ss2)));
title('log10 of mean abs. diff div. by mass')
xlabel('slice')
ylabel('sequence')

subplot(2,1,2);
imagesc(disp_mat');
title('Highlighted selected slices')
xlabel('slice')
ylabel('sequence');





%put the full matrix A together and overwrite t_select with the right hand
%side vector to the equation
I = eye(s);
A = zeros(s*num*(k-s+1),k);
for i = 1:k-s+1
    theta = angles(i:i+s-1);
    Gam = [cos(theta) , sin(theta)];
    Gam = Gam*pinv(Gam)-I;
    for j = 1:num
        t_select(s*(j-1)+1:s*j,i) = ...
            -Gam*t_select(s*(j-1)+1:s*j,i);
        A((i-1)*num*s+s*(j-1)+1:(i-1)*num*s+s*j,i:i+s-1)=Gam;
    end
end

%Shifts solution given by least squares solution
%Careful, if the problem is too big you should use an iterative solver
%instead of the pseudo inverse
shifts = pinv(A)*t_select(:);


%get the shifted stack
for i = 1:k
    stack(:,:,i) = circshift(stack(:,:,i),[round(shifts(i)),0]);
end



        
    

end






















