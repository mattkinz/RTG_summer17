function [slices,init,W,bb,T,v1,v2,chunksizes,numchunks,given]...
    = getDIPS(stack,opt,init)


%Written by: Toby Sanders
%Comp. & Applied Math Dept., Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016



if ~sum(sum(sum(init)))~=0
    if opt.initialsolution
        if ischar(opt.initialfile)
            init=load(opt.initialfile);
        end
        init=struct2cell(init);
        if max(size(init))>1
            fprintf('matlab file for the initial solution has multiple variables\n');
            fprintf('will recompute the initial solution\n');
            init=0;
        else
            init=init{1};
        end
        if size(init,3)~=size(stack,2) && size(init,3)~=opt.endslice-opt.startslice+1
            fprintf('the size of the initial solution doesnt agree with the input size from\n')
            fprintf('startslice and endslice.  Therefore, it cannot be determined\n')
            fprintf('which slices the initial solution are for.  Will compute a new initial solution\n')
            init=0;
        elseif size(init,3)~=opt.endslice-opt.startslice+1
            init = init(:,:,opt.startslice:opt.endslice);
        end
    else
        init=0;
    end
end

oldsizestack = size(stack,2);
stack = stack(:,opt.startslice:opt.endslice,:);




if size(stack,3)~=1
    [numray,slices,angles]=size(stack);
    if angles~=max(size(opt.angles))
        error('the specified angles do not match the size of the stack');
    end
    bb = zeros(numray*angles,slices);
    for i = 1:slices
        bb(:,i)=reshape(stack(:,i,:),numray*angles,1);
    end
else
    slices=size(stack,2);
    bb=stack;
    numray = size(bb,1)/max(size(opt.angles));
    angles = max(size(opt.angles));
end


if ~isfield(opt,'W')
    W = radonmatrix(opt.angles,opt.resolution,numray);
elseif ~sum(sum(opt.W))
    W = radonmatrix(opt.angles,opt.resolution,numray);
elseif (size(opt.W,1)~=numray*angles) || (size(opt.W,2)~=opt.resolution^2)
    fprintf('projection matrix does not match data,\n constructing new matrix\n');
    W = radonmatrix(opt.angles,opt.resolution,numray);
else
    W = opt.W;
end
opt.W='';




[v1,v2,chunksizes]=determinechunks(opt.chunksize,opt.overlap,slices);
numchunks=max(size(chunksizes));



if sum(sum(sum(init)))~=0
    if size(init,3)~=slices && size(init,3)~=oldsizestack
        error('the size of the stack doesnt agree with the initial input solution');
    elseif size(init,3)==oldsizestack
        init = init(:,:,opt.startslice:opt.endslice);
    end
    if size(init,1)*size(init,2)~=size(W,2)
        error('dimension of the initial solution doesnt agree with the specified reconstruction size');
    end
    given=1;
else
    init=zeros(opt.resolution,opt.resolution,slices);
    given=0;
end

if opt.rx~=0
    T = computegaussian(opt.rx,opt.rz,opt.sigmax,opt.sigmaz,slices,opt.chunksize);
else
    T=1;
end


end