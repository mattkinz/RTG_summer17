function forward = forward_proj(rec,angles,numray)

%OUTPUT: forward projections of "rec" at the angles "angles"

%INPUT: 
    %rec - the reconstruction to be forward projected
    %angles - list of the angles to forward project to, in degrees
    %numray - dectector count
    
%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



if size(angles,1)~=1 && size(angles,2)~=1
    W=angles;
    numangle = size(W,1)/numray;
else
    W = radonmatrix(angles,size(rec,1),numray,0);
    numangle = max(size(angles));
end

rec = reshape(rec,size(rec,1)*size(rec,2),size(rec,3));

bb = W*rec;

forward = zeros(numray,size(rec,2),numangle);
for i = 1:size(rec,2)
    forward(:,i,:)=reshape(bb(:,i),numray,1,numangle);
end

end