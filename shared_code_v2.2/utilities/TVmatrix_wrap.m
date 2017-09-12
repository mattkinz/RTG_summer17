function M = TVmatrix_wrap(d1,d2)

%constructs sparse matrix with +1 and -1 to compute 2-D gradient of an
%image of size d1xd2


%builds the TV matrix along the y-dimension, d1
if d1>1, D1 = def_diff_matrix_y(d1,d2); 
else D1=zeros(0,d1*d2); end;

%builds the TV matrix along the x-dimension, d2
if d2>1, D2 = def_diff_matrix_x(d1,d2); 
else D2=zeros(0,d1*d2); end; 

M = [D1;D2];

end


function D1 = def_diff_matrix_y(d1,d2)

y = zeros(2*d1*d2,3);
c=1;
for i = 1:d2
    for j = 1:d1
        if j~=d1
            y(2*(c-1)+1,:)=[c,c,-1];
            y(2*c,:)=[c,c+1,1];
        else
            y(2*(c-1)+1,:)=[c,c,-1];
            y(2*c,:)=[c,c+1-d1,1];
        end
        c=c+1;
    end
end

D1 = sparse(y(:,1),y(:,2),y(:,3),d1*d2,d1*d2);

end


function D2 = def_diff_matrix_x(d1,d2)


c=1;
y = zeros(2*d1*d2,3);


for i = 1:d2-1
    for j = 1:d1
            y(2*(c-1)+1,:)=[c,c,-1];
            y(2*c,:)= [c,c+d1,1];
        c=c+1;
    end
end

for j = 1:d1
   y(2*(c-1)+1,:) = [c,c,-1];
   y(2*c,:) = [c,c-d1*(d2-1),1];
   c=c+1;
end

D2 = sparse(y(:,1),y(:,2),y(:,3),d1*d2,d1*d2);

end




