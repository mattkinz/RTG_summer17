function [U,out,init] = PDIRT(opt,bb,init)


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

opt.grays = [opt.grays(1), mean(opt.grays),opt.grays(end)];

[opt,bb] = checkdirtoptions(opt,bb);

opt.grays(2)='';


[slices,init,W,Wn,bb,c,n,T,v1,v2,chunksizes,numchunks,out,given] = ...
    getDIRT(bb,opt,init);







U = zeros(n,n,slices);





for k = 1:numchunks
    
    %compute an initial solution if needed
    if ~given
        fprintf('\n\nFinding inital SIRT solution for chunk %i of %i\n',k,numchunks);
        init(:,:,v1(k,1):v1(k,2)) = sirtden(bchunk,W,n,60,grays(1));
    end

    %set up chunk stuff
    bchunk = bb(:,v1(k,1):v1(k,2));
    ppp = ceil(chunksizes(k)/2);
    nn = norm(bchunk(:,ppp));
    Uchunk = init(:,:,v1(k,1):v1(k,2));
    
    
    
    %Display the initial stuff
    Ucheck = Uchunk;
    if opt.disppic
        figure(1);
        subplot(2,2,1);imshow((Uchunk(:,:,ppp)-opt.grays(1))/opt.grays(end));
        title('Initial Solution');
        figure(1);
    end
    out.iterative_error(1,k)=norm(W*reshape(Uchunk(:,:,ppp),n^2,1)...
        -bchunk(:,ppp))/nn;
    
    if chunksizes(k)>1
        fprintf('\n\nBeginning refinements for chunk number %i of %i\n',...
            k,numchunks);
        fprintf('This chunk covers slices %i to %i\n\n',v1(k,1),v1(k,2));
    end
    
    
    
    
    for i = 1:opt.total_iter
        

        %compute the smoothing and thresholding 
        %find the refinement pixels and check for convergence
        Uchunk = imfilter(Uchunk,T,'replicate');
        [Uchunk,set,fixed]=thresh_tol_special(Uchunk,opt.thresh,opt.grays);
              
        out.set(i+1,k)=max(size(set{ppp}));
        
        
        
        
        
        %print the iteration information
        if opt.disp
        fprintf('iter=%i, ||Wx-b||/||b||=%g pixels=%i \n',...
            i-1,out.iterative_error(i,k),out.set(i+1,k));
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
        Uchunk = DIRT_update(W,Wn,bchunk,Uchunk,set,fixed,opt.update_type,opt,c,n);

        
    end
    
    
    
    
    %Smooth the overlapping regions together
    for j = v1(k,1):v1(k,2)
        U(:,:,j)=U(:,:,j)+Uchunk(:,:,j-v1(k,1)+1)*v2(j,k);
    end
    out.proj_error(k)=out.iterative_error(end,k);
    
    
    
    
    %print the final iteration info
    if opt.disp
    fprintf('iter=%i, ||Wx-b||/||b||=%g pixels=%i \n',...
        i,out.iterative_error(i+1,k),out.set(i+1,k));
    end
    
    
end
       


%finalize the solution with one last smoothing
out.mean_error_entire = sum(out.proj_error)/numchunks;
out.opt=opt;

U = imfilter(U,T);
U = thresh_tol_special(U,opt.thresh,opt.grays);

end