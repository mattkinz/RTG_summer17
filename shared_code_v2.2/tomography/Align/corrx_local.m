function [stack,v,stackw] = corrx_local(stack,xrange,stackw)


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
startslice=ceil(k/2);
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
for i = startslice-1:-1:1
    T1 = stackw(:,xrange(1)-p:xrange(2)+p,i);
    T2 = stackw(:,xrange(1)-p:xrange(2)+p,i+1);
    X = sum(T1(:,:));
    Y = sum(T2(:,p+1:end-p));
    s=0;
    for j = -p:p
            resid = norm(abs(X(j+p+1:end+j-p)-Y));
        if resid<s || j==-p
            s = resid;
            v(i)=-j;
        end
    end
end



for i = startslice+1:k
    %[T1,T2]=localtrunc(stack(:,xrange(1)-p+1:xrange(2)+p,i-1:i),0);
    T1 = stackw(:,xrange(1)-p:xrange(2)+p,i-1);
    T2 = stackw(:,xrange(1)-p:xrange(2)+p,i);
    X = sum(T2(:,:));
    Y = sum(T1(:,p+1:end-p));
    s=0;
    for j = -p:p
            resid = sum(abs(X(j+p+1:end+j-p)-Y));
        if resid<s || j==-p
            s = resid;
            v(i)=-j;
        end
    end
end


for i = startslice-1:-1:1
   v(i)=sum(v(i:i+1));
   move = v(i);
   stack(:,:,i) = circshift(stack(:,:,i),[0,round(move)]);
   stackw(:,:,i) = circshift(stackw(:,:,i),[0,round(move)]);
end
  

for i = startslice+1:k
   v(i)=sum(v(i-1:i));
   move = v(i);
   stack(:,:,i)=circshift(stack(:,:,i),[0,round(move)]);
   stackw(:,:,i) = circshift(stackw(:,:,i),[0,round(move)]);
end