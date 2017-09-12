function [U,out] = inpaint_3D_joint_batch(bb,S,d,opts)



% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 09/22/2016



if numel(d)~=3
    error('this is 3D inpainting bro!');
end
%p = d(1);q = d(2); r = d(3);

[~,s] = sort(S(:,2));
S(:,1) = S(s,1);S(:,2) = S(s,2);S(:,3) = S(s,3);
bb = bb(s);

chunksize = 50;
overlap = 10;
slices = d(2);
[v1,v2,chunksizes] = determinechunks(chunksize,overlap,slices);
numchunks = numel(chunksizes);
Uchunk = cell(numchunks);
out = cell(numchunks);
d2 = d;


U = zeros(d);
my_timer = zeros(numchunks,1);
for i = 1:numchunks
    tic;
    fprintf('Inpainting chunk %d of %d\n',i,numchunks);
    s = find(v1(i,1)<=S(:,2)  & S(:,2)<=v1(i,2));
    d2(2) = chunksizes(i);
    

    T = sub2ind(d2,S(s,1),S(s,2)-v1(i,1)+1,S(s,3));
    b2 = bb(s);


    if numel(T)~=numel(b2)
        error('number of data points and specified indices dont match');
    end


    A = @(x,mode)subdata_select(x,mode,T,d2);

    [Uchunk{i},out{i}] = HOTV_joint(A,b2,d2,opts);
    
    
    
    my_timer(i) = toc;
    estimated_time = (sum(my_timer)/i*numchunks-sum(my_timer))/60;
    fprintf('\n Estimated remaining time is %g minutes\n',estimated_time);

end



for i = 1:numchunks
    %Uchunk{i}=reshape(Uchunk{i},dim(1),dim(2),chunksizes(i));
    for j = v1(i,1):v1(i,2)
        U(:,j,:)=U(:,j,:)+Uchunk{i}(:,j-v1(i,1)+1,:)*v2(j,i);
    end
end




function x = subdata_select(x,mode,S,d)

switch mode
    case 1
        x = x(S);
    case 2
        y = zeros(d);
        y(S) = x;
        x = y(:);
end