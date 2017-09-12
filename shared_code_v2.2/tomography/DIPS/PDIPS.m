function [U,U2,out] = PDIPS(opts,stack,init)


%Written by: Toby Sanders
%Comp. & Applied Math Dept., Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016

if ~exist('stack','var')
    stack = 0;
end
if ~exist('init','var')
    init=0;
end




%%%Get the options and the set up DIPS
[opts,stack]=check_pdips_options(opts,stack);

[slices,init,W,bb,T,v1,v2,chunksizes,numchunks,given] = ...
    getDIPS(stack,opts,init);
n=opts.resolution;
if ~given
    init=TV3D_auto(bb,W,n,opts);
end




%Initialize the reconstructions
U = zeros(n,n,slices);
out.iter_error = zeros(opts.dips_iter+1,numchunks);
out.set = zeros(opts.dips_iter+1,numchunks);
out.set(:,1)=n^2;

%Main loop that performs the function over all of the chunks
for k = 1:numchunks
        
    
        Uchunk = init(:,:,v1(k,1):v1(k,2));
        bchunk = bb(:,v1(k,1):v1(k,2));
        nn = norm(bchunk,'fro');
        ppp = ceil(chunksizes(k)/2);
        
        
        %get the initial error and print
        if nn~=0
            out.iter_error(1,k) = norm(bchunk...
                -W*reshape(Uchunk,n^2,chunksizes(k)),'fro')/nn;
        end
        fprintf('\n Beginning refinements for chunk %d of %d\n',k,numchunks);


        %loop runs over the prespecified number of updates
    for i = 1:opts.dips_iter
        
        
        %smoothing and segmenting
        Uchunk = imfilter(Uchunk,T);
        [Uchunk,set,fixed]=thresh_tol_special(Uchunk,...
            opts.thresh,opts.grays);
        opts.init=Uchunk;
        out.set(i+1,k)=max(size(set{ppp}));
        fprintf('iter=%i, ||Wx-b||/||b||=%g , pixels=%i\n',...
            i,out.iter_error(i,k),out.set(i+1,k));
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
            s = s + max(size(fixed{j}));
        end
        
        
        %This vector is now set up to include the "fixed" pixels
        %A factor of 10 is used for the fixed scaling
        b3=zeros(size(bb,1)*chunksizes(k)+s,1);
        c=0;
        lambda2 = 2000/opts.poly.mu;
        
        for j = 1:chunksizes(k)
            c2 = c + size(bb,1)+max(size(fixed{j}));
            b3(c+1:c2) = [b2(:,j);Uchunk(fixed{j},j)*lambda2];
            c=c2;
        end
        
        
        
       
        
        %Set up the operater for this refinement and run
        A = @(x,mode)sub_oper(W,x,mode,set,fixed,s,lambda2);
    	Uchunk=poly_ann_3D(A,b3,n,n,chunksizes(k),opts.poly);
        
        
        %get the new error and print
        if nn~=0
            out.iter_error(i+1,k) = norm(bchunk...
                -W*reshape(Uchunk,n^2,chunksizes(k)),'fro')/nn;
        end

    end

    %Put the results of the chunk into the main U
    for j = v1(k,1):v1(k,2)
            U(:,:,j)=U(:,:,j)+Uchunk(:,:,j-v1(k,1)+1)*v2(j,k);
    end
    
end
U2 =U;

if opts.bdry_dirt_iter + opts.region_dirt_iter>0
    fprintf('\n BEGINNING DART ITERATIONS \n');
    [U,out_dirt,~] = DIRT(opts,bb,U);
    out.dirt = out_dirt;
else
    %Smooth once more and segment
    T = computegaussian(opts.rx,opts.rz,2*opts.sigmax,2*opts.sigmaz,slices,opts.chunksize);
   % U = threshold(imfilter(U,T),opts.thresh,opts.grays);
    out.proj_error = norm(W*reshape(U,n^2,slices)-bb);
    out.dirt = 'no DART iterations';
end

end
