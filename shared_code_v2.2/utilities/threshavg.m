function [avg,U] = threshavg(U,threshes,grays)



%THRESHAVG%

%DESCRIPTION:
    %threshavg performs the same function as threshold, but additionally
    %returns the average of all pixels between each of the thresholds

%NOTATION:
    %[avg,X_thresh] = threshavg(X,thresholds,grays)

%INPUTS: 
    %X - any matrix, 1-D, 2-D, 3-D, or N-D, to be thresholded.
    %thresholds - a vector holding the thresholds for rec
    %grays - a vector holding the gray values to be used for the 
        %thresholding.  This variable is optional for this function.

%OUTPUTS: 
    %avg - a vector of the same length as grays.
    %avg(1) is the average of all the pixels less than thresholds(1)
    %avg(i) is the average of all pixels between thresholds(i-1) and
    %thresholds(i)
    %avg(end) is the average of all pixels greater than all of the input
    %thresholds
    %X_thresh - the thresholded reconstruction.
    
    
    
%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


num=max(size(threshes));
if ~exist('grays','var')
    grays = 0:1/(num):1;
end

if max(size(grays))~=num+1
    error('grays should be one larger than threshes');
end

[m,n,k]=size(U);
avga=zeros(k,num+1);
U = reshape(U,m*n,k);
for i = 1:k
    [s,s2]=sort(U(:,i));
    indexold=0;
    for p = 1:num
        indexup=m*n;
        indexdown=indexold+1;
        index=round((indexup+indexdown)/2);
        while indexup-indexdown~=1
            if s(index)>threshes(p)
                indexup=index;
                index=round((index+indexdown)/2);
            else
                indexdown=index;
                index=round((index+indexup)/2);
            end
        end
        U(s2(indexold+1:index),i)=grays(p);
        avga(i,p)=mean(s(indexold+1:index));
        indexold=index;
        if indexold==m*n,break;end;
    end
    avga(i,end)=mean(s(index+1:end));
    U(s2(index+1:end),i)=grays(end);
end
avg=zeros(1,num+1);
for i = 1:num+1
    avg(i)=sum(avga(:,i))/nnz(avga(:,i));
end

U = reshape(U,m,n,k);

end
            
            