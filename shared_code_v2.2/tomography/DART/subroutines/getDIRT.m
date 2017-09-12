function [slices,init,W,d1,bb,c,n,T,v1,v2,chunksizes,numchunks,out,given]...
    = getDIRT(stack,opt,init)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014




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
    W = radonmatrix(opt.angles,opt.resolution,numray,0);
elseif ~sum(sum(opt.W))
    W = radonmatrix(opt.angles,opt.resolution,numray,0);
elseif (size(opt.W,1)~=numray*angles) || (size(opt.W,2)~=opt.resolution^2)
    fprintf('projection matrix does not match data,\n constructing new matrix\n');
    W = radonmatrix(opt.angles,opt.resolution,numray,0);
else
    W = opt.W;
end
opt.W='';




if sum(strcmp('SIRT',{opt.region_update_type,opt.bdry_update_type}))
    %Wn = W';
    %for i = 1:size(Wn,2)
    %    Wn(:,i)=Wn(:,i)/sum(Wn(:,i));
    %end
    %Wn=Wn';

    %c = zeros(size(W,2),1);
    %for i = 1:size(W,2)
    %    c(i)=1/sum(Wn(:,i));
    %end
    d1 = 1./sum(W,2);
    %d1 = repmat(d1,1,slices);
    c = 1./sum(W,1)';
    %c = repmat(c,1,slices);
else
    Wn=[];
    c=[];
end


n = opt.resolution;
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
    init=zeros(n,n,slices);
    given=0;
end



if opt.rx~=0
    T = computegaussian(opt.rx,opt.rz,opt.sigmax,opt.sigmaz,slices,opt.chunksize);
else
    T=1;
end


    
    
out.proj_error = zeros(numchunks,1);
out.set = zeros(opt.total_iter+1,numchunks);
out.set(1,:)=n^2;
out.iterative_error = zeros(opt.total_iter+1,numchunks);

end