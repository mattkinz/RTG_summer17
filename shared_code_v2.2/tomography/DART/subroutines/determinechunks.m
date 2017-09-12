function [v1,v2,chunksizes]= determinechunks(chunksize,overlap,slices)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

%Last update: 02/2016


if slices <= chunksize
    v1=[1,slices];
    v2 = ones(slices);
    chunksizes=slices;
    return;
end

if overlap == 0 
    numchunks = ceil(slices/chunksize);
    v1 = zeros(numchunks,2);
    for i = 1:numchunks-1
        v1(i,1)=(i-1)*chunksize+1;
        v1(i,2)=i*chunksize;
    end
    v1(end,1)=v1(end-1,2)+1;
    v1(end,2)=slices;
    chunksizes = zeros(numchunks,1);
    chunksizes(1:numchunks-1)=chunksize;
    chunksizes(numchunks)=slices-v1(end,1)+1;
    v2 = zeros(slices,numchunks);
    for i = 1:numchunks
        v2(v1(i,1):v1(i,2),i)=1;
    end
    return;
end


if overlap >chunksize || overlap==1 
    overlap=2;
    if chunksize==2
        chunksize=3;
    end
    fprintf('your overlap is poor, adjusting a bit');
end






numchunks = ceil(slices/(chunksize-overlap));

chunksize = round((slices-overlap)/numchunks)+overlap;
v1 = zeros(numchunks,2);
chunksizes = zeros(numchunks,1);
v2 = zeros(slices,numchunks);
for i = 1:numchunks
    v1(i,1)=1+(chunksize-overlap)*(i-1);
    v1(i,2)=(chunksize-overlap)*(i-1) + chunksize;
    if i==numchunks
        v1(i,2)=slices;
    end
    chunksizes(i)=v1(i,2)-v1(i,1)+1;
end

v2(1:v1(2,1),1)=1;
v2(v1(2,1):v1(1,2),1)=1:-1/(overlap-1):0;
for i = 2:numchunks
    if i ~=numchunks
        v2(v1(i-1,2):v1(i+1,1),i)=1;
        v2(v1(i+1,1):v1(i,2),i)=1:-1/(overlap-1):0;
    end
    v2(v1(i,1):v1(i-1,2),i)=0:1/(overlap-1):1;
end
v2(v1(numchunks-1,2):end,end)=1;



end