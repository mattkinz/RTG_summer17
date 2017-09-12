function [W] = radon_alt(angles,dim,numray,scale)


% Written by Toby Sanders @ASU
% Department of Statistical & Mathematical Sciences
% 05/17/2016


% A alternative code for generating the tomography matrix. The design is
% much simpler than the original code (radonmatrix.m), but requires
% significantly more computational work

if ~exist('scale')
    scale=1;
end

if scale>1
    d0 = dim*scale;    
else
    d0 = dim;
end

if scale<1
    d0 = dim*scale;
end



N = max(size(angles));

%distance between parallel rays
H = d0/numray;

%distance of each parallel ray from the center
t = (-d0/2+H/2):H:(d0/2-H/2);


%find the rays that actually intersect with the grid based off of scaling
useable_rays = find(abs(t) < sqrt(2)*dim);
numray_new = size(useable_rays,2);

%image grid point x & y
xpoints = -dim/2+1/2:dim/2-1/2;
ypoints = -dim/2 + 1/2:dim/2-1/2;
[X,Y] = meshgrid(xpoints,ypoints);

accum = [];
for i = 1:N
    a = angles(i);
    abs_a = abs(a);

    %distance that determines if grid point is interested by ray
    limit_distance = 1/2*( cosd(a) + sind(abs_a) );
    
    %compute the distance of each rotated grid point from x-axis
    xprod = X*cosd(a);
    yprod = Y*sind(a);
    Z = xprod+yprod;

    

    useable_rays = find(abs(t) < dim*(abs(sind(a))+abs(cosd(a))));
    numray_new = size(useable_rays,2);
    
    
    %determine which which pixels and rays intersect for this angle
    %set pixi_count to know how many pixels 
    S = cell(numray_new,1);
    pixi_count = zeros(numray_new,1);
    for j = useable_rays
        Zj = Z - t(j);
        [Sy,Sx] = find( (-limit_distance<= Zj) &  (Zj< limit_distance) , 2*dim );
        pixi_count(j) = size(Sy,1);
     
        S{j} = zeros(pixi_count(j),3);
        S{j}(:,1) = dim*(Sx-1) + Sy;
        S{j}(:,2) = abs( Z(S{j}(:,1)) - t(j) );
        
    end
    
    %accumulate all of the ray information into a single vector matrix
    y = zeros(sum(pixi_count),4);
    counter = [0;cumsum(pixi_count)];
    for j = useable_rays
        y(counter(j) + 1:counter(j+1),1) = (i-1)*numray + j;
        y(counter(j) + 1:counter(j+1),2) = S{j}(:,1);
        y(counter(j) + 1:counter(j+1),3) = S{j}(:,2);
    end
    
    
    %Find and store the lengths of all intersections, based on the 2 cases
    %The S1 case is when the ray passes all the way through the cell (in
    %the top and out the bottom, for example)
    %The S2 case is otherwise
    S1 = 1:sum(pixi_count);
    if abs_a<=45
        S2 = find(y(:,3) > 1/2*( cosd(a) - sind(abs_a) ) );
        S1(S2) = '';
        y(S1,4) = abs(secd(a));%secd(a);
    else
        S2 = find(y(:,3) > 1/2*( cosd(90-abs_a) - sind(abs(90-abs_a)) ) );
        S1(S2) = '';
        y(S1,4) = abs(cscd(a));% cscd(abs_a);
    end
    y(S2,4) = (limit_distance - y(S2,3))/(abs(cosd(a)*sind(a)));

    %accumulate all of the information into the matrix-vector accum
    accum = [accum;y(:,[1,2,4])];    
end

W = sparse(accum(:,1),accum(:,2),accum(:,3),numray*N,dim^2);