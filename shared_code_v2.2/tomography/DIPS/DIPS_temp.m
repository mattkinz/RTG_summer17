function [U,U2,out] = DIPS(opts,stack,init)


%Written by: Toby Sanders
%Comp. & Applied Math Dept., Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016

%OUTPUTS:
    %U - the 3-D final segmented dips reconstruction.
    %U2 - the result from DIPS without final segmentation or DART
    %iterations
    %out - data about the reconstruction

    
    
if ~exist('stack','var')
    stack = 0;
end
if ~exist('init','var')
    init=0;
end




%%% Get the options and set up DIPS
[opts,stack]=check_dips_options(opts,stack);

[slices,init,W,bb,T,v1,v2,chunksizes,numchunks,given] = ...
    getDIPS(stack,opts,init);
n=opts.resolution;
if ~given
    init= HOTV3D_tomo(bb,W,n,opts);
end



%Initialize variables
U = zeros(n,n,slices);
out.iter_error = zeros(opts.dips_iter+1,numchunks);
out.set = zeros(opts.dips_iter+1,numchunks);
out.set(:,1)=n^2;
out.r_tol = zeros(opts.dips_iter+1,numchunks);
out.run_time = zeros(numchunks,1);

%Main loop that performs the function over all of the chunks
for k = 1:numchunks
        
    tic;
        Uchunk = init(:,:,v1(k,1):v1(k,2));
        bchunk = bb(:,v1(k,1):v1(k,2));
        nn = norm(bchunk,'fro');
        ppp = ceil(chunksizes(k)/2);
        
        
        %get the initial error and print
        if nn~=0
            out.iter_error(1,k) = norm(bchunk...
                -W*reshape(Uchunk,n^2,chunksizes(k)),'fro')/nn;
        end
        fprintf('\n Beginning refinements for block %d of %d\n',k,numchunks);

        
        
        c_tol=opts.t_tol;
        out.r_tol(1,k) = c_tol;
        %loop runs over the prespecified number of updates
    for i = 1:opts.dips_iter
        
        
        %smooth and segment current solution
        Uchunk = imfilter(Uchunk,T,'replicate');
        [Uchunk,set,fixed]=thresh_tol(Uchunk,...
            opts.grays,c_tol);
        opts.poly.init=Uchunk;
        out.set(i+1,k)=max(size(set{ppp}));
        fprintf('iter=%i, ||Wx-b||/||b||=%g , pixels=%i, r_tol=%g\n',...
            i,out.iter_error(i,k),out.set(i+1,k),c_tol);
        %display the refinement region
       %{     
        if opts.disppic
            temp = zeros(n);
            temp(set{ppp})=1;
            figure(2);imshow(temp);figure(2);
        end
        %}
        
        %setting up the vector for updating
        b2 = zeros(size(bb,1),chunksizes(k));
        Uchunk = reshape(Uchunk,n^2,chunksizes(k));
        
        s=0;
        for j = 1:chunksizes(k)
            b2(:,j) = bchunk(:,j) - W(:,fixed{j})*Uchunk(fixed{j},j);
            s = s + numel(fixed{j});
        end
        
        
        
        %This vector is now set up to include the "fixed" pixels
        %A factor of 10 is used for the fixed scaling
        %b3=zeros(size(bb,1)*chunksizes(k)+s,1);
        %c=0;
        %lambda2 = 2000/opts.poly.mu;
        
        %for j = 1:chunksizes(k)
        %    c2 = c + size(bb,1)+max(size(fixed{j}));
        %    b3(c+1:c2) = [b2(:,j);Uchunk(fixed{j},j)*lambda2];
        %    c=c2;
        %end
        b3 = [];
        for j = 1:chunksizes(k)
            b3 = [b3;Uchunk(fixed{j},j)];
        end
        col_fixed = zeros(s,1);
        c = 0;
        for j = 1:chunksizes(k)
            c2 = numel(fixed{j});
            col_fixed(c+1:c2) = fixed{j};
            c = c2;
        end

        
        
       
        
        %Set up the operater for this refinement and run
        if i>1
            opts.poly.sigma = pa_out.sigma;
            opts.poly.delta = pa_out.delta;          
        end
        A = @(x,mode)sub_oper_temp(W,x,mode,set);
    	[Uchunk,pa_out]=DIPSTV(A,b2,b3,col_fixed,[n,n,chunksizes(k)],opts.poly);
        
        
        %get the new error and print
        if nn~=0
            out.iter_error(i+1,k) = norm(bchunk...
                -W*reshape(Uchunk,n^2,chunksizes(k)),'fro')/nn;
        end
        
        
        
        if abs(out.set(i+1,k)-out.set(i,k))/out.set(i+1,k)<opts.t_epsilon
            c_tol = c_tol + opts.t_delta;
            if c_tol>max(diff(opts.grays))/2
                fprintf('r_tol is now too high, moving on\n\n');
                break;                
            end        
        end
        out.r_tol(i+1,k) = c_tol;
        %c_tol = c_tol+(opts.f_tol-opts.t_tol)/(opts.dips_iter-1);
    end

    %Put the results of the chunk into the main U
    for j = v1(k,1):v1(k,2)
            U(:,:,j)=U(:,:,j)+Uchunk(:,:,j-v1(k,1)+1)*v2(j,k);
    end
    
    %print remaining time
    out.run_time(k)=toc;
    estimated_time = (sum(out.run_time)/k*numchunks-sum(out.run_time))/60;
    fprintf('\n Estimated remaining time is %g minutes\n',estimated_time);
end
U2 =U;

if opts.bdry_dirt_iter + opts.region_dirt_iter>0
    fprintf('\n BEGINNING DART ITERATIONS \n');
    [U,out_dirt,~] = DIRT(opts,bb,U);
    out.dirt = out_dirt;
else
    %Smooth once more and segment
    U = threshold(imfilter(U,T,'replicate'),opts.thresh,opts.grays);
    out.proj_error = norm(W*reshape(U,n^2,slices)-bb);
    out.dirt = 'no DART iterations';
end

end
