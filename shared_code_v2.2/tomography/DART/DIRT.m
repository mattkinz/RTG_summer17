function [U,out,init] = DIRT(opt,bb,init)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

%Last updated on 05/2016.

%OUTPUTS:
    %U - the 3-D DART reconstruction.  The slices change with the 3rd
        %dimension, i.e. dartrec(:,:,i) holds the ith slice of the
        %reconstruction
    %out - data about the reconstruction
    %init - the SIRT reconstruction or input reconstruction used as
        %the initial solution to compute "dartrec".


%INPUTS:
%Inputs that must be user defined:
    %opt.angles - the projection angles listed in order, in degrees
    %opt.thresh - a vector holding the thresholds used for segmenting 
    %opt.grays - a vector holding the gray values used in segmentation.
    %opt.data_name - name of the file in which the aligned tilt series is 
        %saved.  Note the tilt axis should be horizontal and centered.
    %opt.data_type - type of file of opt.filename, set to either 'mrc' or 
        %'mat'

%Inputs recommended to be user defined, but default values are otherwise
%used:
    %opt.bdry_dirt_iter - number of DART iterations with boundary cleaning
    %opt.region_dirt_iter - number of region refinements (DIP-LS)
    %opt.t_tol - radius of interval used for the partial segmentation
        %function with DIPS.  Set a small tolerance for slower more
        %careful convergence and a high tolerance for fast convergence
    %opt.t_delta - increase in opt.t_tol whenever new pixels are not being
        %classified
    %opt.t_epsilon - this epsilon is used to determine if opt.t_tol should 
        %be increased by opt.t_delta.  It is the case if the number of 
        %classified pixels has only changed relatively by less that
        %epsilon.
    %opt.bdry_update_type - indicates the solver used for bdry updates.
        %Set to 'SIRT' (default) or 'CGLS'.
    %opt.region_update_type - indicates the solver used for region updates.
        %Set to 'SIRT' (default) or 'CGLS'. 
    %opt.rand_pixel_probability - probability to a free pixel is set to
        %free.  Default is 0.
    %opt.resolution - the number of pixels used for reconstruction.
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
    %opt.inner_iter - number of SIRT iterations in each DART iteration.
        %Default value is 10.
    %opt.W - if the projection matrix is precomputed, set opt.W to be
        %the projection matrix.
    %opt.disp - prints information about the reconstrution at each
        %iteration.  Set to "true" for the display and "false" 
        %otherwise.  Default is "true".
    %opt.disppic - displays various images about the reconstruction at
        %each iteration. Set to "true" for this display and "false"
        %otherwise.  Default is "true".
    %opt.convergence_criteria - if set to "true", DART will iterate until
        %the reconstruction converges.  Set to "false" to simply use
        %opt.outeriter instead.  Default is "false".  Typically
        %convergence is seen after just a few iterations.
    %opt.convergence_tol - tolerance for opt.convergencecriteria.
        %Default is 10^(-3).
    %opt.max_iter - maximum number of iterations if
        %opt.convergencecriteria is used.  Default is 50.
    %opt.rx - radius of the gaussian in the slice (x,y) dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.rz - radius of the gaussian in the slice z dimension used
        %for smoothing.  Set to 1 for no smoothing.  Default is 3.
    %opt.sigmax - sigma value used for the gaussian smoothing in the
        %(x,y) dimension.  Default is 1.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.sigmaz - sigma value used for the gaussian smoothing in the
        %z dimension.  Default is 1.5.  Increase for more smoothing
        %and decrease for less smoothing.
    %opt.chunksize - size of the chunks that DART will run on.  Default
        %is 50.
    %opt.overlap - once the chunks are computed, they are merged
        %together using a partition of unity.  The overlap is the overlap
        %of this merging.  Default is 10.






if ~exist('bb','var')
    bb=false;
end

if ~exist('init','var')
    init=false;
end

[opt,bb] = checkdirtoptions(opt,bb);

[slices,init,W,Wn,bb,c,n,T,v1,v2,chunksizes,numchunks,out,given] = ...
    getDIRT(bb,opt,init);



%%USER OPTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%thresh=opt.thresh;
%grays=opt.grays;
opt.bad_pixels=false;
region_dirt_iter = opt.region_dirt_iter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



U = zeros(n,n,slices);




out.run_time = zeros(numchunks,1);
for k = 1:numchunks
    
    tic;
    %compute an initial solution if needed
    bchunk = bb(:,v1(k,1):v1(k,2));
    if ~given
        fprintf('\n\nFinding inital SIRT solution for block %i of %i\n',k,numchunks);
        init(:,:,v1(k,1):v1(k,2)) = SIRT(bchunk,W,n,30,opt.grays(1));
    end

    %set up chunk stuff
    
    ppp = ceil(chunksizes(k)/2);
    nn = norm(bchunk(:,ppp));
    Uchunk = init(:,:,v1(k,1):v1(k,2));
    
    
    
    %Display the initial stuff
    Ucheck = Uchunk;
    if opt.disppic
        figure(1);
        subplot(2,2,1);imagesc(Uchunk(:,:,ppp),[opt.grays(1),opt.grays(end)]);
        title('Initial Solution');
        figure(1);
    end
    out.iterative_error(1,k)=norm(W*reshape(Uchunk(:,:,ppp),n^2,1)...
        -bchunk(:,ppp))/nn;
    
    if chunksizes(k)>1
        fprintf('\n\nBeginning refinements for block number %i of %i\n',...
            k,numchunks);
        fprintf('This block covers slices %i to %i\n\n',v1(k,1),v1(k,2));
    end
    
    
    
    %Reset the values which were possibly changed in previous chunk
    opt.region_dirt_iter = region_dirt_iter;
    c_tol=opt.t_tol;
    
    for i = 1:opt.total_iter
        

        %compute the smoothing and thresholding 
        %find the refinement pixels and check for convergence
        Uchunk = imfilter(Uchunk,T,'replicate');
        [Uchunk,Ucheck,set,fixed,update_type,type,convergence] = ...
            DIRT_thresh(Uchunk,opt,i,c_tol,Ucheck);
        if convergence
            break;
        end        
        out.set(i+1,k)=max(size(set{ppp}));
        
        
        
        
        
        %print the iteration information
        if opt.disp
        fprintf('iter=%i, ||Wx-b||/||b||=%g pixels=%i,type=%s \n',...
            i-1,out.iterative_error(i,k),out.set(i+1,k),type);
        end
        
        
        
        
        
        
        %display images
        if opt.disppic
            DIRT_display(Uchunk(:,:,ppp),opt,out,i,k);
        end
        
        
        
        %Compute the error after thresholding
        Uchunk = reshape(Uchunk,n^2,chunksizes(k));
        ee = norm(bchunk(:,ppp)-W*Uchunk(:,ppp));
        if nn~=0
            out.iterative_error(i+1,k) = ee/nn;
        end
        
        
        
        

        
        %Perform the updating
        Uchunk = DIRT_update(W,Wn,bchunk,Uchunk,set,fixed,update_type,opt,c,n);

        
        %Check for appropriate r_tol and convergence of the solution
        if i<opt.region_dirt_iter
            if abs(out.set(i+1,k)-out.set(i,k))/out.set(i+1,k)<opt.t_epsilon
                c_tol = c_tol + opt.t_delta;
                if c_tol>min(diff(opt.grays))/2
                    opt.region_dirt_iter=i;
                    fprintf('r_tol is now too high, switching to boundary\n\n');
                end        
            end
        end
    end
    
    
    
    
    %Smooth the overlapping regions together
    for j = v1(k,1):v1(k,2)
        U(:,:,j)=U(:,:,j)+Uchunk(:,:,j-v1(k,1)+1)*v2(j,k);
    end
    out.proj_error(k)=out.iterative_error(end,k);
    
    
    
    
    %print the final iteration info
    if opt.disp
    fprintf('iter=%i, ||Wx-b||/||b||=%g pixels=%i,type=%s \n',...
        i,out.iterative_error(i+1,k),out.set(i+1,k),type);
    end
    
    
    %calculate time remaining
    out.run_time(k)=toc;
    estimated_time = (sum(out.run_time)/k*numchunks-sum(out.run_time))/60;
    fprintf('\n Estimated remaining time is %g minutes\n',estimated_time);
end



%finalize the solution with one last large smoothing
out.mean_error_entire = sum(out.proj_error)/numchunks;
out.opt=opt;

%{
if opt.rx~=0
    T = computegaussian(opt.rx*2,opt.rz*2,opt.sigmax*2,...
        opt.sigmaz*2,slices,opt.chunksize);
else
    T=1;
end
%}

U = imfilter(U,T,'replicate');
U = threshold(U,opt.thresh,opt.grays);

end