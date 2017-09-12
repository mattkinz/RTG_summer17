function [stack,v,stackw] = corrx_global(stack,xrange,stackw)

%CORRX_GLOBAL and CORRX_LOCAL

%DESCRIPTION: 
    %these functions performs alignment in the "x-direction," i.e. alignment
    %shifts parallel to the tilt axis.  The alignment criteria is based on
    %the conservation of mass in each slice.  In other words, the mass in a
    %cross section or slice perpendicular to the tilt axis should say in 
    %that slice for all tilts.
    
%NOTATION:
    %[stacknew,s,stackwnew] = corrx_global(stack,xrange,stackw);
    %[stacknew,s,stackwnew] = corrx_local(stack,xrange,stackw);
    
%INPUTS:
    %stack - the tilt series, with the tilt axis horizontal
    %xrange - a vector with 2 entries, which tells which cross sections of
        %the stack to start and end at for the alignment computations
    %stackw - an optional input windowing of the tilt series.  If the user
        %inputs "stackw", then this tilt series is for the computations to
        %determine the shifts for alignment of "stack."  If this input is
        %skipped, then a windowed tilt series "stackw" is generated within
        %the algorithm and then used for alignment.
        
%OUTPUTS:
    %stacknew - the tilt series, after the alignment
    %s - the horizontal shifts that were applied
    %stackwnew - the windowed tilt series, after the alignment.
    
    
%What's the difference?
    %Both of these algorithms perform alignment based on the idea that the 
    %mass along cross sections perpendicular to the tilt axis remains
    %constant (idealy).  "corrx_global" computes the global average of each
    %cross section(the average mass of each cross section of each 
    %projection) , and aligns each projection based on this global average.
    % "corrx_local" computes the mass along each cross section of a single
    % projection and aligns the adjacent projection based on the single
    % previous (or prior) projection, hence the name local.
    
%NOTES:
    %CORRX_GLOBAL is more susceptible to local error but should produce a
    %good global alignment, i.e. the beginning and ending projections
    %should be well aligned, but there may be small errors between adjacent
    %projections.  The opposite holds for CORRX_LOCAL.  CORRX_GLOBAL is
    %recommended.
    
    

%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014
[~,n,k]=size(stack);
p=20;
if xrange ==0
    xrange=[round(n/5),round(4/5*n)];
else
    xrange=[round(xrange(1))+p,round(xrange(2))-p];
end
if ~exist('stackw','var')
    stackw = weighting2(stack,'super');
end



v = zeros(k,1);
X = zeros(k,xrange(2)-xrange(1)+2*p+1);
for i = 1:k
    X(i,:)=sum(stackw(:,xrange(1)-p:xrange(2)+p,i));
end
mm = mean(X(:,p+1:end-p));
figure(2);subplot(1,2,1);imagesc(X);
title('sums of the projections');figure(2);

for i = 1:k
    s =0;
    for j = -p:p
    resid = norm(mm-X(i,p+j+1:end-p+j));
        if resid<s || j==-p
            s = resid;
            v(i)=-j;
        end
    end
end

for i = 1:k
   stack(:,:,i) = circshift(stack(:,:,i),[0,v(i)]);
   stackw(:,:,i) = circshift(stackw(:,:,i),[0,v(i)]);
   X(i,:)=circshift(X(i,:),[0,v(i)]);
end

figure(2);subplot(1,2,2);imagesc(X);
title('sums of the projections after alignment');figure(2);