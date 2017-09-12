
function [U,out] = HOTV3D_tomo(stack,angles,recsize,opts)


% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


%PA3D_tomo

%DESCRIPTION:
    %this function performs 3 dimensional HOTV minimization
    %reconstructions of tomography data.  It makes use of the HOTV3D code
    %originally called the TVAL3 code
    %written by Chengbo Li, alumni of Rice University.
    
%NOTATION:
    %rec = PA3D_tomo(stack,angles,recsize,opts);
    
%INPUTS:
    %stack - the tilt series, where it is assumed the tilt axis is 
    %horizontal and located at the middle of the stack
    %angles - a vector holding the projection angles of the stack, in order,
    %in degrees
    %recsize - the dimension of the reconstruction
    %   opts: structure containing input parameters, 
    %       see function check_HOTV_opts.m for these



chunksize = 20;
overlap = 3;
slices = size(stack,2);
dim = recsize;
if max(size(dim))==1
    dim(1) = recsize;
    dim(2) = recsize;
end

U = zeros(dim(1),dim(2),slices);
[v1,v2,chunksizes]=determinechunks(chunksize,overlap,slices);
numchunks=max(size(chunksizes));
Uchunk = cell(numchunks);

if size(angles,1)==1 || size(angles,2)==1
    if size(stack,3)==1
        numray = size(stack,1)/numel(angles);
    elseif size(stack,3)~=numel(angles)
        error('Number of input angles doesnt match the number of projection images');
    else
        numray = size(stack,1);
    end
    W = radonmatrix(angles,dim,numray,0);
else
    W = angles;
end

my_timer = zeros(numchunks,1);
out = cell(numchunks,1);

if isfield(opts,'reweighted_TV')
    if opts.reweighted_TV
        w = reshape(opts.coef,dim(1),dim(2),slices,3);
        init = opts.init;
    end
else
    opts.reweighted_TV = false;
end


for i = 1:numchunks
    tic;
    fprintf('Reconstruction of chunk %d of %d\n',i,numchunks);
    num_slice = v1(i,2)-v1(i,1)+1;
    
    bb = zeros(size(stack,1)*size(stack,3)*num_slice,1);
    for j = 1:num_slice
        bb((j-1)*size(stack,1)*size(stack,3)+1:j*size(stack,1)*size(stack,3))=...
            col(stack(:,v1(i,1)+j-1,:));
    end
    if opts.reweighted_TV
        opts.coef = reshape(w(:,:,v1(i,1):v1(i,2),:),dim(1)*dim(2)*num_slice,3);
        opts.init = init(:,:,v1(i,1):v1(i,2));
    end
    A = @(x,mode)mm_operator(x,mode,W,num_slice,dim);
    
    [Uchunk{i},out{i}]=HOTV3D(A,bb,[dim(1),dim(2),chunksizes(i)],opts);
    
    out{i}.chunk_range = [v1(i,1),v1(i,2)];
    my_timer(i)=toc;
    estimated_time = (sum(my_timer)/i*numchunks-sum(my_timer))/60;
    fprintf('\n Estimated remaining time is %g minutes\n',estimated_time);
end




for i = 1:numchunks
    Uchunk{i}=reshape(Uchunk{i},dim(1),dim(2),chunksizes(i));
    for j = v1(i,1):v1(i,2)
        U(:,:,j)=U(:,:,j)+Uchunk{i}(:,:,j-v1(i,1)+1)*v2(j,i);
    end
end




function y = mm_operator(x,mode,W,slices,res)

switch mode
    case 1
        x = reshape(x,res(1)*res(2),slices);
        y = col(W*x);
    case 2
        y = col((reshape(x,size(W,1),slices)'*W)');
end
