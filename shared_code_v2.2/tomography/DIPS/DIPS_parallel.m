function [U,U2,out] = DIPS_parallel(opts,stack,init)

if ~exist('stack','var')
    stack = 0;
end
if ~exist('init','var')
    init=0;
end




%%%Get the options and the set up DIPS
[opts,stack]=check_dips_options(opts,stack);

[slices,init,W,bb,T,v1,v2,chunksizes,numchunks,given] = ...
    getDIPS(stack,opts,init);
n=opts.resolution;
if ~given
    init=TV3D_auto(bb,W,n,opts);
end




%Initialize the reconstructions
out.iter_error = zeros(opts.dips_iter+1,numchunks);
out.set = zeros(opts.dips_iter+1,numchunks);
out.set(:,1)=n^2;
out.r_tol = zeros(opts.dips_iter,numchunks);


%Main loop that performs the function over all of the chunks

%my_timer = zeros(numchunks,1);
Uchunk = cell(numchunks,1);
bchunk = cell(numchunks,1);
S = cell(numchunks,1);
outvar = cell(numchunks,1);

for k = 1:numchunks
    Uchunk{k} = init(:,:,v1(k,1):v1(k,2));
    bchunk{k} = bb(:,v1(k,1):v1(k,2));
    S{k} = opts;
    outvar{k} = out;
end

clear init;

parfor k = 1:numchunks
        
   % tic;
       % Uchunk{k} = init(:,:,v1(k,1):v1(k,2));
        
        nn = norm(bchunk{k},'fro');
        ppp = ceil(chunksizes(k)/2);
        
        
        %get the initial error and print
        if nn~=0
            outvar{k}.iter_error(1,k) = norm(bchunk{k}...
                -W*reshape(Uchunk{k},n^2,chunksizes(k)),'fro')/nn;
        end
        fprintf('\n Beginning refinements for chunk %d of %d\n',k,numchunks);

        
        
        c_tol=S{k}.t_tol;
        outvar{k}.r_tol(1,k) = c_tol;
        %loop runs over the prespecified number of updates
    for i = 1:S{k}.dips_iter
        
        
        %smoothing and segmenting
        Uchunk{k} = imfilter(Uchunk{k},T,'replicate');
        [Uchunk{k},set,fixed]=thresh_tol(Uchunk{k},...
            S{k}.grays,c_tol);
        S{k}.init=Uchunk{k};
        outvar{k}.set(i+1,k)=max(size(set{ppp}));
        
        fprintf('iter=%i, ||Wx-b||/||b||=%g , pixels=%i, r_tol=%g\n',...
            i,outvar{k}.iter_error(i,k),outvar{k}.set(i+1,k),c_tol);
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
        Uchunk{k} = reshape(Uchunk{k},n^2,chunksizes(k));
        
        s=0;
        for j = 1:chunksizes(k)
            b2(:,j) = bchunk{k}(:,j) - W(:,fixed{j})*Uchunk{k}(fixed{j},j);
            s = s + max(size(fixed{j}));
        end
        
        
        %This vector is now set up to include the "fixed" pixels
        %A factor of 10 is used for the fixed scaling
        b3=zeros(size(bb,1)*chunksizes(k)+s,1);
        c=0;
        lambda2 = 2000/S{k}.mu;
        
        for j = 1:chunksizes(k)
            c2 = c + size(bb,1)+max(size(fixed{j}));
            b3(c+1:c2) = [b2(:,j);Uchunk{k}(fixed{j},j)*lambda2];
            c=c2;
        end
        
        
        
       
        
        %Set up the operater for this refinement and run
        A = @(x,mode)sub_oper(W,x,mode,set,fixed,s,lambda2);
    	Uchunk{k}=TVAL3D(A,b3,n,n,chunksizes(k),S{k});
        
        
        %get the new error and print
        if nn~=0
            outvar{k}.iter_error(i+1,k) = norm(bchunk{k}...
                -W*reshape(Uchunk{k},n^2,chunksizes(k)),'fro')/nn;
        end
        
        
        
        if abs(outvar{k}.set(i+1,k)-outvar{k}.set(i,k))/outvar{k}.set(i+1,k)<S{k}.t_epsilon
            c_tol = c_tol + S{k}.t_delta;
            if c_tol>max(diff(S{k}.grays))/2
                fprintf('r_tol is now too high, moving on\n\n');
                break;                
            end        
        end
        outvar{k}.r_tol(i+1,k) = c_tol;
        %c_tol = c_tol+(opts.f_tol-opts.t_tol)/(opts.dips_iter-1);
    end

    
    Uchunk{k}=single(Uchunk{k});
    %print remaining time
    %my_timer(k)=toc;
    %estimated_time = (sum(my_timer)/k*numchunks-sum(my_timer))/60;
    %fprintf('\n Estimated remaining time is %g minutes\n',estimated_time);
end



%Put the results of the chunk into the main U
U = zeros(n,n,slices,'single');

for k = 1:numchunks
    for j = v1(k,1):v1(k,2)
            U(:,:,j)=U(:,:,j)+Uchunk{k}(:,:,j-v1(k,1)+1)*v2(j,k);
    end
    
    for i = 1:S{k}.dips_iter+1
       out.set(i,k) = outvar{k}.set(i,k);
       out.iter_error(i,k) = outvar{k}.iter_error(i,k);  
       out.r_tol(i,k) = outvar{k}.r_tol(i,k);
    end
end
clear Uchunk;

U2 =U;

if opts.bdry_dirt_iter + opts.region_dirt_iter>0
    fprintf('\n BEGINNING DART ITERATIONS \n');
    U = double(U);
    [U,out_dirt,~] = DIRT_temp(opts,bb,U);
    out.dirt = out_dirt;
else
    %Smooth once more and segment
   % T = computegaussian(opts.rx,opts.rz,2*opts.sigmax,2*opts.sigmaz,slices,opts.chunksize);
    U = threshold(imfilter(U,T,'replicate'),opts.thresh,opts.grays);
    out.proj_error = norm(W*reshape(U,n^2,slices)-bb);
    out.dirt = 'no DART iterations';
end

end
