function [] = COMdetermine(stack,angles)

%INPUTS 
    %stack - tilt series with horizontal tilt axis in the middle of the
    %stack
    %angles - vector containing the projection6 angles of each projection
    %of the tilt series, in degrees

%This function may be used to check the quality of the center of mass
%alignment.  The center of mass of any slice should follow a
%circular path around the tilt axis.  So for each 2-D slice, the center
%of masses of all the 1-D projections of this slice are computed and plotted 
%together, and then a best fit "circular path" approximation is computed 
% to all of the computed center of masses.  If the alignment and stack are 
%nice, then this approximating curve should fit the computed curve nicely.
%Whether or not this fitting is accurate should only be taken into
%consideration for the slices which contain significant mass.  To see all
%of the plot, just run this function After a few seconds a plot should show 
%up for the first slice, just press any keystroke to continue for all of 
%the other slices


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

    [k,which] = max(size(angles));
    if which ==2,angles=angles';end
    angles = angles*pi/180;
    Phi = zeros(k,2);
    Phi(:,1)=cos(angles);
    Phi(:,2)=-sin(angles);
    
    [m,n,k2]=size(stack);
    
    if k2~=k
        error('stack size doesnt match angles');
    end
    
    
    %stack = weighting2(stack,'super');
    t=zeros(k,n);
    w = (1:m)';
    
    for i = 1:n
        for j = 1:k
            t(j,i)=sum(w.*stack(:,i,j))/sum(stack(:,i,j));
        end
    end
    t = t-(m+1)/2;
    s = zeros(k,1);
    for i =1:n
        c(:,i)=Phi\t(:,i);
        r = Phi*c(:,i);
        s = s+(reshape(sum(stack(:,i,:)),k,1).*((r-t(:,i)).^2));
    end
        
   %{     
    for i = 1:k
        r = Phi(i,:)*c;
        d(i,:) = r-t(i,:);
        shift(i)=mean(d(i,:));
        st(i)=std(d(i,:));
        for j =1:n
            if abs(d(i,j)-shift(i))>st(i)
                d(i,j)=0;
            end
        end
        shift2(i)=sum(d(i,:))/nnz(d(i,:));
    end
    %}

    for i =1:n
        c = Phi\t(:,i);
        s = Phi*c;
        c
        plot(angles*180/pi,t(:,i),'kx',angles*180/pi,s,'k-','linewidth',2)
        title(['slice ',num2str(i)]);
        %legend('Calculated COM','Least squares fit','Location','NorthWest');
        ylabel('COM location');xlabel('Projection angle, \theta')
        pause;
    end


    
end
    
    









