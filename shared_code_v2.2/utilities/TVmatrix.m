function M = TVmatrix(d1,d2)

%constructs sparse matrix with +1 and -1 to compute 2-D gradient of an
%image of size d1xd2

y = zeros(2*(2*d1*d2-d1-d2),3);
c=1;
for i = 1:d2
    for j = 1:d1-1
        y(2*(c-1)+1,:)=[c,c+(i-1),-1];
        y(2*c,:)=[c,c+i,1];
        c=c+1;
    end
end

c2=1;
for i = 1:d2-1
    for j = 1:d1
        y(2*(c-1)+1,:)=[c,c2,-1];
        y(2*c,:)= [c,c2+d1,1];
        c=c+1;
        c2=c2+1;
    end
end

M = sparse(y(:,1),y(:,2),y(:,3),(d1-1)*d2+(d2-1)*d1,d1*d2);

end




