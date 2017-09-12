function U = threshold(U,threshes,grays)


%DESCRIPTION:
    %function for thresholding

%NOTATION:
    %X_thresh = threshold(X,thresholds,grays);

%INPUTS:
    %X - any matrix, 1-D, 2-D, 3-D, or N-D, to be thresholded.
    %thresholds - a vector holding the thresholds for rec
    %grays - a vector holding the gray values to be used for the 
        %thresholding

%OUTPUT:
    %X - the thresholded matrix

%Notes:  each entry in X is determined to be between 2 of the input
%thresholds or less than or greater than all of the thresholds.  This entry
%is then assigned to the appropriate gray value.  The thresholds should be
%in increasing value, and the gray values may be set to whatever the user
%likes, as long as the number of gray values is one greater than the number
%of thresholds


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


num=max(size(threshes));
if max(size(grays))~=num+1
    error('grays should be one larger than threshes');
end
       
[m,n,k]=size(U);
U = reshape(U,m*n,k);

if num>1
    a = 9/8*threshes(1)-1/8*threshes(end);
    b = (threshes(end)-a)/.9;
    threshes = (threshes-a)/b;
    U = (U-a)/b;
else
    U = U-threshes+1/2;
    threshes=1/2;
end
Ut = cell(num,1);

for i = 1:num
    Ut{i}=im2bw(U,threshes(i));
end
U = zeros(m*n,k)+grays(1);
for i =1:num
    U = U + double(Ut{i})*(grays(i+1)-grays(i));
end

U = reshape(U,m,n,k);

end
            
            