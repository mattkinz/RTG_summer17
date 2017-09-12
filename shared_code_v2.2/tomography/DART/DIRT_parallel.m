function [U,outs,init] = DIRT_parallel(opt,bb,init)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


%INPUTS:
%Inputs that must be user defined:
    %opt.angles - the projection angles listed in order, in degrees
    %opt.thresh - a vector holding the thresholds used for segmenting 
    %opt.grays - a vector holding the gray values used in segmentation.
    %opt.filename - name of the file in which the aligned tilt series is 
        %saved.  Note the tilt axis should be horizontal and centered.
    %opt.filetype - type of file of opt.filename, set to either 'mrc' or 
        %'mat'

%Inputs recommended to be user defined, but default values are otherwise
%used:
    %opt.bdry_dirt_iter - number of DART iterations, i.e. number of iterations of
        %segmenting and back projecting.  Default is 10.
    %opt.region_dirt_iter - number of DART iterations, using region refinements
    %opt.recsize - the number of pixels used for reconstruction.
        %Reconstruction slices will be opt.recsize by opt.recsize. Default
        %value is the detector count
    %opt.startslice - which slice the reconstruction begins on.  Default is 1. 
    %opt.endslice - which slice the reconstruction ends on.  Default is the
    %end of the stack.
    %opt.initialsolution - Set to "true" an initial solution has been
        %computed and will be used, otherwise set to "false".
    %opt.initialfile - name of the matlab file where the reconstruction is
        %saved
        
%More advanced inputs, and default values are again used these are not
%specified:
    %opt.rescale - after reconstruction, the slices are rescaled to be
        %the number of detector counts, so that when the reconstruction 
        %is viewed in 3-D there's is correct scaling in the 3
        %dimensions.  Set to "true" for the rescaling or "false" is no
        %rescaling is desired.  Default is "true"
    %opt.inneriter - number of SIRT iterations in each DART iteration.
        %Default value is 10.
    %opt.W - if the projection matrix is precomputed, set opt.W to be
        %the projection matrix.
    %opt.disp - prints information about the reconstrution at each
        %iteration.  Set to "true" for the display and "false" 
        %otherwise.  Default is "true".
    %opt.disppic - displays various images about the reconstruction at
        %each iteration. Set to "true" for this display and "false"
        %otherwise.  Default is "true".
    %opt.convergence - if set to "true", DART will iterate until
        %the reconstruction converges.  Set to "false" to simply use
        %opt.outeriter instead.  Default is "false".  Typically
        %convergence is seen after just a few iterations.
    %opt.conv_tol - tolerance for opt.convergencecriteria.
        %Default is 10^(-6).
    %opt.maxiter - maximum number of iterations if
        %opt.convergencecriteria is used.  Default is 50.
    %opt.rx - radius of the gaussian in the slice (x,y) dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.rx - radius of the gaussian in the slice z dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.sigmax - sigma value used for the gaussian smoothing in the
        %(x,y) dimension.  Default is 1.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.sigmaz - sigma value used for the gaussian smoothing in the
        %z dimension.  Default is 1.5.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.chunksize - size of the chunks that DART will run on.  Default
        %is 50.
    %opt.overlap - one the chunks are computed, they are merged
        %together using a partition of unity.  The overlap is the overlap
        %of this merging.  Default is 10.





if ~exist('bb','var')
    bb=false;
end

if ~exist('init','var')
    init=false;
end

[opt,bb] = checkdirtoptions(opt,bb);

[slices,init,W,Wn,bb,c,n,T,v1,v2,chunksizes,numchunks,~,given] = ...
    getDIRT(bb,opt,init);



%%USER OPTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%thresh=opt.thresh;
%grays=opt.grays;
opt.bad_pixels=false;
opt.convergence_criteria=false;
%region_dirt_iter = opt.region_dirt_iter;
total_iter = opt.total_iter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



U = zeros(n,n,slices);


%compute an initial solution if needed
    if ~given
        fprintf('Finding the initial solution \n\n');
        init = sirtauto(bb,W,n,50,grays(1));
    end

bchunk = cell(numchunks,1);
Uchunk = cell(numchunks,1);
outs = cell(numchunks,1);

for k = 1:numchunks
    bchunk{k} = bb(:,v1(k,1):v1(k,2));
    Uchunk{k} = init(:,:,v1(k,1):v1(k,2));
    outs{k} = zeros(total_iter+1,1);
end
c_tol=opt.t_tol;

parfor k = 1:numchunks

    %set up chunk stuff
    %bchunk = bb(:,v1(k,1):v1(k,2));
    ppp = ceil(chunksizes(k)/2);
    nn = norm(bchunk{k}(:,ppp));
    %Uchunk = init(:,:,v1(k,1):v1(k,2));
    
    
    
    %Display the initial stuff
    %{
    Ucheck = Uchunk{k};
    if opt.disppic
        figure(1);
        subplot(2,2,1);imshow((Uchunk(:,:,ppp)-opt.grays(1))/opt.grays(end));
        title('Initial Solution');
        figure(1);
    end
    %}
   % out.iterative_error(1,k)=norm(W*reshape(Uchunk{k}(:,:,ppp),n^2,1)...
    %    -bchunk(:,ppp))/nn;
    outs{k}(1) = norm(W*reshape(Uchunk{k}(:,:,ppp),n^2,1)-bchunk{k}(:,ppp))/nn;
    if chunksizes(k)>1
        fprintf('\n\nBeginning refinements for chunk number %i of %i\n',...
            k,numchunks);
        fprintf('This chunk covers slices %i to %i\n\n',v1(k,1),v1(k,2));
    end
    
    
    
    %Reset the values which were possibly changed in previous chunk
    %opt.region_dirt_iter = region_dirt_iter;
   % 
    
    for i = 1:total_iter
        

        %compute the smoothing and thresholding 
        %find the refinement pixels and check for convergence
        Uchunk{k} = imfilter(Uchunk{k},T,'replicate');
        [Uchunk{k},~,set,fixed,update_type,type,convergence] = ...
            DIRT_thresh(Uchunk{k},opt,i,c_tol,[]);
        if convergence
            break;
        end        
       %out.set(i+1,k)=max(size(set{ppp}));
        
        
        
        
        
        %print the iteration information
       % if opt.disp
        %fprintf('iter=%i, ||Wx-b||/||b||=%g pixels=%i,type=%s \n',...
        %    i-1,out.iterative_error(i,k),out.set(i+1,k),type);
        %end
        
        
        
        
        
        
        %display images
       % if opt.disppic
        %    DIRT_display(Uchunk{k}(:,:,ppp),opt,out,i,k);
        %end
        
        
        
        %Compute the error after thresholding
        Uchunk{k} = reshape(Uchunk{k},n^2,chunksizes(k));
        ee = norm(bchunk{k}(:,ppp)-W*Uchunk{k}(:,ppp));
        if nn~=0
            outs{k}(i+1) = ee/nn;
        end
        
        
        
        

        
        %Perform the updating
        Uchunk{k} = DIRT_update(W,Wn,bchunk{k},Uchunk{k},set,fixed,update_type,opt,c,n);

        
        %Check for appropriate r_tol and convergence of the solution
       %{
        if i<opt.region_dirt_iter
            if abs(out.set(i+1,k)-out.set(i,k))/out.set(i+1,k)<.005
                c_tol = c_tol + opt.t_delta;
                if c_tol>min(diff(opt.grays))/2
                    c_tol=c_tol/2;
                    fprintf('r_tol is now too high, switching to boundary\n\n');
                end        
            end
        end
        %}
        %c_tol = c_tol + (opt.f_tol-opt.t_tol)/(opt.region_dirt_iter-1);
        %if mod(i,3)==0
        %    opt.bad_pixels=false;
        %else
        %    opt.bad_pixels=false;
        %end
    end
    
    
    
    %print the final iteration info
    fprintf('\n Chunk %i of %i complete,\n',k ,numchunks);
    fprintf('||Wx0-b||/||b||=%g \n ||Wx1 - b||/||b||=%g \n,',outs{k}(1),outs{k}(2));
    fprintf('||Wx%i -b||/||b||=%g \n',i,outs{k}(end));
    
end

%Smooth the overlapping regions together
for k = 1:numchunks
    for j = v1(k,1):v1(k,2)
        U(:,:,j)=U(:,:,j)+Uchunk{k}(:,:,j-v1(k,1)+1)*v2(j,k);
    end
end

%finalize the solution with one last large smoothing

if opt.rx~=0
    T = computegaussian(opt.rx*2,opt.rz*2,opt.sigmax*2,...
        opt.sigmaz*2,slices,opt.chunksize);
else
    T=1;
end

U = imfilter(U,T,'replicate');
U = threshold(U,opt.thresh,opt.grays);

end