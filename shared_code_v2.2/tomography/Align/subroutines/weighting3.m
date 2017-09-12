function [stackw ] = weighting3(stack,accuracy )

if ~exist('accuracy','var')
    accuracy = 'medium';
elseif ~ischar(accuracy)
    accuracy = 'medium';
end

[m,n,k]=size(stack);
stackw = zeros(m,n,k);

chunksize=10;
overlap = 0;

[v1,v2,chunksizes]=determinechunks(chunksize,overlap,n);
numchunks=max(size(chunksizes))

for i =1:numchunks
    i
    s=weighting2(stack(:,v1(i,1):v1(i,2),:),accuracy);
    for j =v1(i,1):v1(i,2)
        stackw(:,j,:)=stack(:,j,:)+s(:,j-v1(i,1)+1,:)*v2(j,i);
    end
end

end

