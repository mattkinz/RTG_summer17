function W = radonmatrix(angles,dim,numray,weights,scale,model)


%DESCRIPTION: 
    %This function generates a projection or Radon matrix, that can
    %generate projections from input projection geometry.
    
%NOTATION:
    %W = radonmatrix(angles,dim,numray,scale,model);

%INPUTS:
    %angles - a vector containing the projection angles in degrees
    %dim - the dimension of the 2-D functions (or images) that this matrix
        %will be able to project, i.e. this matrix will project images that
        %will pixel dimensions dimxdim.
    %numray - number of projection measurements per angle.
    %scale - a real positive number which scales the projection process to 
        %only a fraction of the image.  For example, if scale = 1/2 and 
        %d =512, then the detection start at 129 and end at 256+128.  This
        %input is optional and the default is no scaling.
    %model - determines whether to use the line or strip model for the
        %projections.  Input 'line' for the line model and 'strip' for the
        %strip model.  Default is 'line'.  I think the 'strip' model has
        %some bugs as well

%OUTPUT: Projection matrix W


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%Mathematics Department, Arizona State University
%10/23/2015



if min(angles)<-89.5 || max(angles)>360
    error('Angles should be contained within the interval (-89.5,360)');
end

if max(size(weights))~=numray*max(size(angles))
    error('weights do not match specified numray and angles');
end



if ~exist('scale','var')
    scale=1;
    remove_pixel_flag=0;
    dimf=max(dim);
elseif scale>1
   dimf=dim;
   dim = round(dim*scale);
   scale=1;
   remove_pixel_flag=1;
else
    remove_pixel_flag=0;
    dimf=max(dim);
end




if ~exist('model','var')
    model = 'line';
elseif ~ischar(model)
    model = 'line';
elseif ~sum(strcmp(model,{'line','strip'}))
    model = 'line';
end




fprintf('\nCalculating projection weights for matrix\n')
if strcmp(model,'line')
    W = radon_line(angles,dim,numray,weights,scale,dimf,remove_pixel_flag);
else
    W = radon_strip(angles,dim,numray,scale);
end





function W = radon_line(angles, d, numray,weights,scale,dimf,remove_pixel_flag)



if max(size(d))==1
    d1=d;
    d2=d;
else
    d1=d(1);
    d2=d(2);
end
d=max(d1,d2);

N = max(size(angles));

if scale>1 || scale <0
    scale=0;
end


if scale~=0
    d0=scale*d;
else
    d0=d;
end

y = zeros(2*d*N*numray,3,'double');
x = zeros(2*d,2,'single');

for i = 1:N

    a = angles(i);
    beta_flag = 0;
    ninety_flag=0;
    gamma_flag=0;
    
    
        
        
        if a>90 && a<270
            a = a-180;
            beta_flag=1;
        elseif a>270.5
            a = a-360;
        elseif abs(270-a)<.5
            a = a-180;
            beta_flag=1;
        end
        
        if abs(90-abs(a))<.5
            ninety_flag=1;
            a = 90-a;
        end
        
        if a<0
            a = abs(a)*pi/180;
            gamma_flag=1;
        else
            a = a*pi/180;
        end

        
        
        H=d0/numray;
        tann =tan(a);
        secc=sec(a);
        cscc=csc(a);
        cott=cot(a);
        s=H*secc;
        L = sqrt((d/2+cos(pi/2+a)*d/2)^2/sin(pi/2+a)^2+d^2/4);
        p1 = cot(asin(d/2/L))*d/2;
        p2 = cot(pi/2-a/2)*d/2;
        p=sin(a)*(p1+p2);

        b=0;
        p=p+s/2+(d-d0)/(2*sin(pi/2-a));
    for j = 1:numray
        %this is to adjust for the next ray and check if we crossed the
        %corner
        if p<d && b==0
            entry=1; 
            h = p-floor(p);
            index = d*floor(p)+1;
        elseif p<d && b==1;
            entry=0;
            h=p-floor(p);
            index = d*d-d+ceil(p)+floor(1-(ceil(p)-p));
        else
            b = 1;
            entry=0;
            p = cott*(p-d);
            h=p-floor(p);
            s = H*cscc;
            index = d*d-d+ceil(p)+floor(1-(ceil(p)-p));
        end
        
        x(:,1)=0;
        x(:,2)=1;
        for k = 1:2*d
            
            if entry ==1    
                if tann<h
                    entry =1;
                    x(2*d-k+1,1)=secc;
                    h = h-tann;
                    indexn = index+1;
                elseif tann>h
                    entry = 0;
                    x(2*d-k+1,1)=h*cscc;
                    h = h*cott;
                    indexn = index-d;
                else
                    entry = 1;
                    h=1;
                    x(2*d-k+1,1)=secc;
                    indexn=index-d+1;
                end
                x(2*d-k+1,2)=index;
            else
                if cott<1-h
                    entry =0;
                    h=h+cott;
                    x(2*d-k+1,1)=cscc;
                    indexn = index-d;
                elseif cott>1-h
                    entry=1;
                    x(2*d-k+1,1)=secc*(1-h);
                    h = 1-(1-h)*tann;
                    indexn = index+1;
                else
                    entry=1;
                    h=1;
                    x(2*d-k+1,1)=cscc;
                    indexn = index-d+1;
                end
                x(2*d-k+1,2)=index;
            end
           
              if floor(index/d)==index/d, break;end;
              if indexn<1, break;end
              if indexn>d*d, break;end
              index=indexn;
        end
            %reweighting
            x(:,1) = x(:,1)*weights((i-1)*numray+j);
            
            
            i1 = (i-1)*numray*2*d+(j-1)*2*d+1;
            i2 = (i-1)*numray*2*d+j*2*d;
            y(i1:i2,3)=x(:,1);
            
            
            
            
        if ~beta_flag
            y(i1:i2,1)=(i-1)*numray+j;
        else
            y(i1:i2,1) = (i-1)*numray+numray-j+1;
        end
        
        
        
        if gamma_flag
            y(i1:i2,2) = d-mod(x(:,2)-1,d)+floor((x(:,2)-1)/d)*d;            
        else
            y(i1:i2,2)=x(:,2);
        end
        
        if ninety_flag
            y(i1:i2,2) = floor((y(i1:i2,2)-1)/d)+1+mod(y(i1:i2,2)-1,d)*d;
        end
     p=p+s;
    end  


end


if remove_pixel_flag
   fprintf('\nRemoving pixels out of bounds\n')
   y = remove_pixels(y,d,dimf); 
end

fprintf('\nBuilding sparse matrix\n');
W = sparse(y(:,1),y(:,2),y(:,3),numray*N,dimf*dimf)*numray/d0;

if ~remove_pixel_flag
    l1(1)=(dimf-d1)/2+1;l1(2)=(dimf-d1)/2+d1;
    l2(1)=(dimf-d2)/2+1;l2(2)=(dimf-d2)/2+d2;
    pixels = get_pixels(l1,l2,dimf);
    W = W(:,pixels);
end



function y = remove_pixels(y,dim,dimf)
    j0 = floor((dim-dimf)/2)+1;
    j1 = floor((dim-dimf)/2)+dimf;
    
    j = uint16(mod(y(:,2),dim));
    k = uint16(floor(y(:,2)/dim)+1);
    j = j.*uint16(j>j0 & j<j1);
    k = k.*uint16(k>j0 & k<j1);
    
    ind = find((j~=0)&(k~=0));
    j = j(ind);j = j- floor((dim-dimf)/2);
    k = k(ind);k = k - floor((dim-dimf)/2);
    y = y(ind,:);
    
    y(:,2) = single(j)+single((k-1))*dimf;
    
   
    
function pixels = get_pixels(level1,level2,mainsize)

m = level1(2)-level1(1)+1;
n = level2(2)-level2(1)+1;

pixels=zeros(m*n,1);
for i = 1:n
    loc1 = (level2(1)+i-2)*mainsize+level1(1);
    pixels((i-1)*m+1:i*m) = loc1:loc1+m-1;
end



function W = radon_strip(angles, d, numray,scale)


if max(size(d))==1
    d1=d;
    d2=d;
end
d=max(d1,d2);

N = max(size(angles));
d0=scale*d2;



y = zeros(2*d*N*numray*5,3,'double');
x = zeros(2*d*5,2);

for i = 1:N


    a = angles(i);
        if a<0
            a = abs(a)*pi/180;
            gamma=0;
        else
            a = a*pi/180;
            gamma=1;
        end
        
        H=d0/numray;
        tann =tan(a);
        secc=sec(a);
        cscc=csc(a);
        cott =cot(a);
        s=H*secc;
        L = sqrt((d/2+cos(pi/2+a)*d/2)^2/sin(pi/2+a)^2+d^2/4);
        p1 = cot(asin(d/2/L))*d/2;
        p2 = cot(pi/2-a/2)*d/2;
        p=sin(a)*(p1+p2);
        b0=0;
        p=p+s/2+(d-d0)/(2*sin(pi/2-a));
        p0=p;
        s0=s;
    for j = 1:numray
        x(:,1)=0;
        x(:,2)=1;
        if p0<d2 && b0==0

            elseif p<d && b0==1;

            else
                b0 = 1;
                p0 = cott*(p0-d);
                s0 = H*cscc;
        end
        b=b0;
        p=p0;
        if a==0
            l1 = p-d0/numray/2;
            l2 = p+d0/numray/2;
            r1 = floor(l1);
            r2 = ceil(l2);
            if r2-r1==1
                x(1:d,1)=1;
                x(1:d,2)=(r1)*d+1:(r1+1)*d;
            else
                for k = r1:r2-1
                    k0=k-r1;
                    if k==r1
                        x(k0*d+1:(k0+1)*d,1)=1-(l1-r1);
                        x(k0*d+1:(k0+1)*d,2)=k*d+1:(k+1)*d;
                    elseif k==r2-1
                        if r2 ~=l2;
                            x(k0*d+1:(k0+1)*d,1)=1-(r2-l2);
                            x(k0*d+1:(k0+1)*d,2)=(k)*d+1:(k+1)*d;
                        end
                    else
                        x(k0*d+1:(k0+1)*d,1)=1;
                        x(k0*d+1:(k0+1)*d,2)=k*d+1:(k+1)*d;
                    end
                end
            end
        else
            for rays=1:5
                if b==0
                    p = p0+(rays-3)*secc/4*d0/numray;
                else
                    p = p0+(rays-3)*cscc/4*d0/numray;
                end
            %this is to adjust for the next ray and check if we crossed the
            %corner

                if p<0
                    p=-tann*p;
                    entry=1; 
                    h = p-floor(p);
                    index = d*floor(p)+1;
                elseif p<d && b==0
                    entry=1; 
                    h = p-floor(p);
                    index = d*floor(p)+1;
                elseif p<d && b==1;
                    entry=0;
                    h=p-floor(p);
                    index = d*d-d+ceil(p)+floor(1-(ceil(p)-p));
                else
                    b = 1;
                    entry=0;
                    p = cott*(p-d);
                    h=p-floor(p);
                    index = d*d-d+ceil(p)+floor(1-(ceil(p)-p));
                end
                for k = 1:2*d

                    
                    if entry ==1    
                        if tann<h
                            entry =1;
                            x(2*d-k+1+(rays-1)*2*d,1)=secc;
                            h = h-tann;
                            indexn = index+1;
                        elseif tann>h
                            entry = 0;
                            x(2*d-k+1+(rays-1)*2*d,1)=h*cscc;
                            h = h*cott;
                            indexn = index-d;
                        else
                            entry = 1;
                            h=1;
                            x(2*d-k+1+(rays-1)*2*d,1)=secc;
                            indexn=index-d+1;
                        end
                        x(2*d-k+1+(rays-1)*2*d,2)=index;
                    else
                        if cott<1-h
                            entry =0;
                            h=h+cott;
                            x(2*d-k+1+(rays-1)*2*d,1)=cscc;
                            indexn = index-d;
                        elseif cott>1-h
                            entry=1;
                            x(2*d-k+1+(rays-1)*2*d,1)=secc*(1-h);
                            h = 1-(1-h)*tann;
                            indexn = index+1;
                        else
                            entry=1;
                            h=1;
                            x(2*d-k+1+(rays-1)*2*d,1)=cscc;
                            indexn = index-d+1;
                        end
                        x(2*d-k+1+(rays-1)*2*d,2)=index;

                    end
                      if floor(index/d)==index/d, break;end;
                      if indexn<1, break;end
                      if indexn>d*d, break;end
                      index=indexn;
                end
            end
        end
            i1 = (i-1)*numray*2*d*5+(j-1)*2*d*5+1;
            i2 = (i-1)*numray*2*d*5+j*2*d*5;
            if a==0
                y(i1:i2,3)=x(:,1)/sum(x(:,1))*d;
            else
                y(i1:i2,3)=x(:,1)/sum(x(:,1))*sum(x(4*d+1:6*d,1));
            end
        if gamma
            y(i1:i2,2)=x(:,2);
            y(i1:i2,1)=(i-1)*numray+j;
        else 
            y(i1:i2,1)=(i-1)*numray+numray-j+1;
               y(i1:i2,2)=max(1,(d/2-floor((x(:,2)-1/2)/d))*d...
                   +d*d/2+(mod(x(:,2),d)-d)...
                   +floor((d-mod(x(:,2),d))/d)*d);               
        end
     p0=p0+s0;
    end  


end

W = sparse(y(:,1),y(:,2),y(:,3),numray*N,d*d)*numray/d0;
l1(1)=(d-d1)/2+1;l1(2)=(d-d1)/2+d1;
l2(1)=(d-d2)/2+1;l2(2)=(d-d2)/2+d2;
pixels = get_pixels(l1,l2,d);
W = W(:,pixels);

    