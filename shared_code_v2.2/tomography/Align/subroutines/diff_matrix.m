function D = diff_matrix(m,n,k)

y = zeros(2*(m-1)*n*k,1);
for j = 1:n
    for i = 1:k
        indexing1 =   2*(m-1)*(i-1);
        indexing2 =   (m)*(i-1);
        indexing3 = (m-1)*(i-1);
        
        y(indexing1+1:indexing1+m/2,1)=indexing3+1:indexing3+m/2;
        y(indexing1+m/2+1:indexing1+m-1,1)=indexing3+m/2+1:indexing3+m-1;
        y(indexing1+m:indexing1+3*m/2-1,1)=indexing3+1:indexing3+m/2;
        y(indexing1+3*m/2:indexing1+2*m-2,1)=indexing3+m/2+1:indexing3+m-1;
        
        y(indexing1+1:indexing1+m-1,2)=indexing2+1:indexing2+m-1;
        y(indexing1+m:indexing1+2*m-2,2)=indexing2+2:indexing2+m;
        
        y(indexing1+1:indexing1+m/2,3)=1;
        y(indexing1+m/2+1:indexing1+m-1,3)=-1;
        y(indexing1+m:indexing1+3*m/2-1,3)=-1;
        y(indexing1+3*m/2:indexing1+2*m-2,3)=1;
    end
end

D = sparse(y(:,1),y(:,2),y(:,3));

        
    