function U = sirtauto(stack,angles,recsize,iter,minc,maxc)


%DESCRIPTION:
    %this function run the SIRT algorithm for tomographic reconstruction.
    
%NOTATION:
    % U = sirtauto(stack,angles,recsize,iterations,minc,maxc);

%INPUTS:    
    % stack - the tilt series, where it is assumed the tilt axis is horizontal 
    %   and located at the middle of the stack
    % angles - a vector holding the projection angles of the stack, in order,
    %   in degrees
    % recsize - the dimension of the reconstruction
    % iter - the number of SIRT iterations
    % minc - minimum density constraint, e.g. U>=0.
    % maxc - maximum density constraint, e.g. U<=1.

    % default values are used for recsize, iter, and minc if they are not
    % specified, therefore one may simply input "sirtden(stack,angles)."
    % the recsize will be set to the detector count, i.e. size(stack,1)
    % 50 iterations is default, and no density constraint is used if it is not
    % specified

%OUTPUT: 
    %U - the reconstruction from the input tilt series and other input 
        %parameters

%What is the difference between SIRTDEN and SIRTAUTO?
%It is useful to perform the back projection process over a chunk of slices
%instead of each individual slice to increase stability in the
%reconstructions.  Therefore the "sirtden" function performs the
%backprojection process over the entire input stack.  However, for more
%accuracy in the reconstruction, it is better to do this over smaller
%chunks, like 50 slices instead of 1000, and then put all of the
%chunk reconstructions together in a nice way.  "sirtauto" does this, so
%sirtauto is convenient whenever one is ready to perform the entire SIRT
%reconstruction of all the slices.

% Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

% Assign default values to unspecified parameters
if (nargin < 3 || isempty(recsize)), recsize = 512; end
if (nargin < 4 || isempty(iter)), iter = 50; end
if (nargin < 5 || isempty(minc)), minc=[]; end
if (nargin < 6 || isempty(maxc)), maxc=[]; end



[W,~,~]=getSIRT(angles,recsize,stack);
slices = size(stack,2);
dim = recsize;
U = zeros(dim,dim,slices);
[v1,v2,chunksizes]=determinechunks(50,5,slices);
numchunks=max(size(chunksizes));
Uchunk = cell(numchunks);


mytimer = zeros(numchunks,1);

for i = 1:numchunks
    tic;
    fprintf('Reconstruction of block %d of %d\n',i,numchunks);
    bchunk = stack(:,v1(i,1):v1(i,2),:);
    

    Uchunk{i}=SIRT(bchunk,W,recsize,iter,minc,maxc);
    mytimer(i) = toc;
    remaining_time = (sum(mytimer)*numchunks/i-sum(mytimer))/60;
    fprintf('estimated remaining time is %g minutes\n',remaining_time);
    %figure(1);imagesc(Uchunk{i}(:,:,ceil(size(bchunk,2)/2)));
    %title('Slice from previous chunk');figure(1);
end

for i = 1:numchunks
    Uchunk{i}=reshape(Uchunk{i},dim,dim,chunksizes(i));
    for j = v1(i,1):v1(i,2)
        U(:,:,j)=U(:,:,j)+Uchunk{i}(:,:,j-v1(i,1)+1)*v2(j,i);
    end
end





end