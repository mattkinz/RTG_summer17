function [W,Q] = radonmatrix_phaseCorr(angles,dim,numray,...
    elevation,Wx,maxWr,minF)


%DESCRIPTION: 
    %This function generates a projection or Radon matrix, that can
    %generate projections from the input projection geometry.
    
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
        %input is optional and the default is no scaling.  If scale is
        %bigger than 1 then it is assumed some rays pass outside the pixel
        %grid.  Set scale to 1 or 0 for no scaling.
%VARIABLES USED FOR PHASE CORRECTION
    %elev - elevation of the plane relative to the scene
    %x_mat - dim x dim matrix in which each entry is the spacial
        %x-coordinate of the scene
    %y_mat - dim x dim matrix in which each entry is the spacial
        %y-coordinate of the scene
    %z_mat - dim x dim matrix in which each entry is the spacial
        %z-coordinate of the scene.  In many cases it's all zeros.
    %minF - minimum frequency from the band data



%OUTPUT: Projection matrix W, with an appropriate phase correction
    %Identifying matrix Q, which tells us which rays passed through the
    %pixelated grid based on the scale variable.  If scale<=1, then all
    %rays pass through the grid and so Q is all ones of size numray x N.
    %This can be used for deleting empty rows of W.


% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 05/17/2016



if min(angles)<-89.5 || max(angles)>360
    error('Angles should be contained within the interval (-89.5,360)');
end







%fprintf('\nCalculating projection weights for matrix\n')

    [W,Q] = radon_line_ph(angles,dim,numray,...
       elevation,Wx,maxWr,minF);






function [W,Q] = radon_line_ph(angles,d,numray,...
    elevation,Wx,maxWr,minF)



if max(size(d))==1
    d1=d;
    d2=d;
else
    d1=d(1);
    d2=d(2);
end
d=max(d1,d2);

N = max(size(angles));

scale = maxWr/Wx/cosd(elevation(1));


%define overall scene range
d0=scale*d;


%distance between parallel rays
H = d0/numray;

%distance of each parallel ray from the center
t = (-d0/2+H/2):H:(d0/2-H/2);


%initialize variables
y = cell(numray,N);
x = zeros(2*d,2,'single');
Q = zeros(numray,N,'int8');
counter = zeros(numray,N,'uint64');


pgrs = waitbar(0,sprintf('Computing matrix values (angle %i of %i)',1,N));
for i = 1:N

    a = angles(i);
    
    %reset flags
    beta_flag = 0;
    ninety_flag=0;
    gamma_flag=0;
    
    
        
        %all of these flags are used to make use of the symmetries in the
        %problem.  The main problem is solved for the case 0<=a<90, and the
        %other cases are found by simply changing the indexing at the end
        if a>=90.5 && a<270
            a = a-180;
            beta_flag=1;
        elseif a>=270.5
            a = a-360;
        elseif abs(270-a)<.5
            a = a-180;
            beta_flag=1;
        end
        
        if abs(90-abs(a))<.5
            ninety_flag=1;
            a = a - 90;
        end
        
        if a<0
            a = abs(a)*pi/180;
            gamma_flag=1;
        else
            a = a*pi/180;
        end

        
        %determine beams with in image mesh
        max_dist = d/2*(sin(a)+cos(a));
        useable_rays = find( abs(t) < max_dist );
        Q(useable_rays,i)=1;
        
        %set constants
        tann =tan(a);
        secc=sec(a);
        cscc=csc(a);
        cott=cot(a);
        
        %moving distance along scene edge
        s=H*secc;
        
        %starting position p
        p = d/2*(tann+1)+t(useable_rays(1))*secc;
        
        
        %b is just an indicator
        b=0;
        
        
    %calcuate weights for each useable ray
    for j = useable_rays
        
        %variables used in this loop:
        %entry - 1 indicates that the entry of the line into the pixel is
        %from the top.  0 indicates that the entry is from the side.
        %p - distance the ray starts on the pixel grid
        %index - index that the ray is passing through and being computed
        %h - indicates where along the pixel index that the ray enters
        %x - stores the values for each ray before being put into y
        
        
        
        
        %adjust for the next ray and check if we crossed the corner
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
        
        %find all of the values for the current ray
        %There are at most 2*d pixels intersected
        for k = 1:2*d
            
            %determine how the ray passes through the pixel, and reset
            %variables accordingly
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
           
            %determine if the ray is at the end of the grid
              if floor(index/d)==index/d, break;end;
              if indexn<1, break;end
              if indexn>d*d, break;end
              index=indexn;
        end
        %store the weights into y
        counter(j,i)=k;
        y{j,i} = zeros(k,3,'single');    
        y{j,i}(:,3) = x(2*d-k+1:end,1);
            
            
        %Store the indexing for the matrix.  Check the flags for the
        %symmetries being used
        if ~beta_flag
            y{j,i}(:,1)=(i-1)*numray+j;
        else
            y{j,i}(:,1) = (i-1)*numray+numray-j+1;
        end
        
        if gamma_flag
            y{j,i}(:,2) = d-mod(x(2*d-k+1:end,2)-1,d)+floor((x(2*d-k+1:end,2)-1)/d)*d;            
        else
            y{j,i}(:,2) =x(2*d-k+1:end,2);
        end
        
        if ninety_flag
            y{j,i}(:,2) = floor((y{j,i}(:,2)-1)/d)+1+mod(y{j,i}(:,2)-1,d)*d;
        end
        
        
     %set new distance along grid for the next ray
     p=p+s;
    end  
%fprintf('Angle %i of %i completed\n',i,N);
waitbar(i/2/(N+1),pgrs,sprintf('Computing matrix values (angle %i of %i)',i,N));
end


% Define speed of light (m/s)
c = 299792458;


%accum will hold the final numbers
numvals = sum(sum(double(counter)));
accum = zeros(numvals,3);


%index keepers
ray_counter = [zeros(1,N);cumsum(counter)];
angle_counter = [0,cumsum(ray_counter(end,:))];



%apply phase correction
%fprintf('\nApplying phase correction\n');

vec = linspace(-Wx/2,Wx/2,d);
[x_mat,y_mat] = meshgrid(vec,vec);
z_mat = zeros(d);

if max(size(elevation))==1
    elevation = elevation*ones(N,1);
end

if max(size(minF))==1
    minF = minF*ones(N,1);
end

for i = 1:N
       
        temp_guy = zeros(ray_counter(end,i),3);
        useables = find(Q(:,i))';
        for j = useables
            temp_guy(ray_counter(j,i)+1:ray_counter(j+1,i),:) = y{j,i};
        end
        a = angles(i);
    

        dR = x_mat * cosd(elevation(i)) * cosd(a) + ...
            y_mat * cosd(elevation(i)) * sind(a) + ...
            z_mat * sind(elevation(i));


        phCorr = exp(-1i*4*pi*minF(i)/c*dR(:));
        temp_guy(:,3) = temp_guy(:,3).*phCorr(temp_guy(:,2));
        accum(angle_counter(i)+1:angle_counter(i+1),:) = temp_guy;
        %fprintf('Angle %i of %i completed\n',i,N);
        waitbar((i+N)/2/(N+1),pgrs,sprintf('Applying phase correction (angle %i of %i)',i,N));    
end
clear y;
waitbar((i+N)/2/(N+1),pgrs,'Building sparse matrix');
%y = reshape(y,2*d*N*numray,3);

%delete empty values in y to clear up memory
%s = find(y(:,3));
%y = y(s,:);


%fprintf('\nBuilding sparse matrix\n');
%W = sparse(y(:,1),y(:,2),y(:,3),numray*N,d*d)*numray/d0;
W = sparse(accum(:,1),accum(:,2),accum(:,3),numray*N,d*d);%*numray/d0;
waitbar(1,pgrs);
close(pgrs);
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




    