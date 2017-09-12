function [U,levels,init,outs] = graysearch(gray1,gray2,gray3,opt,stack)

%This function applies DART to various ranges of gray levels,
%which should be specified in opt.gray1,opt.gray2,opt.gray3
%See the file, dirt_rec_test for the opt, options
%The testing is done over just the middle slice of the stack

%Inputs:the options (see dirt_rec_test.m)
%the initial reconstruction doesn't need to be specified

%Outputs: U holds the reconstructions for the various gray values
%The jth row of the matrix levels 
%tells what grays values are used to reconstruct U(:,:,j)


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014




opt.grays=zeros(1,3);
opt.grays(1)=gray1(1);
opt.grays(2)=gray2(1);
opt.grays(3)=gray3(1);
a = max(size(gray1));
b = max(size(gray2));
c = max(size(gray3));
if ~exist('stack','var')
    stack=false;
end

[opt,stack]=checkdirtoptions(opt,stack);
slices = opt.startslice:opt.endslice;

if ~isfield(opt,'W')
    opt.W = radonmatrix(opt.angles,opt.resolution,size(stack,1),0);
elseif sum(sum(opt.W))==0
    opt.W = radonmatrix(opt.angles,opt.resolution,size(stack,1),0);
end

U = zeros(opt.resolution,opt.resolution,a*b*c);

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
    init = sirtden(stack(:,slices,:),opt.W,opt.resolution,30,gray1(1));
end

outs = cell(a*b*c,1);
levels=zeros(a*b*c,5);
levels(:,1)=1:a*b*c;
for i0 = 1:a
    opt.grays(1)=gray1(i0);
    for i1 = 1:b
        opt.grays(2)=gray2(i1);
        for i2 = 1:c
            opt.grays(3)=gray3(i2);
                fprintf('\n__________________________\n');
                fprintf('gray pair =[%g,%g,%g]\n',opt.grays(1),...
                    opt.grays(2),opt.grays(3));
                num=(i0-1)*b*c+(i1-1)*c+i2;
                [T,out]=...
                    DIRT(opt,stack,init);
                U(:,:,num)=T(:,:,ceil(size(T,3)/2));
                outs{num}=out;
                levels(num,2)=gray1(i0);
                levels(num,3)=gray2(i1);
                levels(num,4)=gray3(i2);
                levels(num,5)=outs{num}.proj_error;
        end
    end
end

init = init(:,:,ceil(size(T,3)/2));

plot(1:a*b*c,levels(:,5),'b-o');
title('errors for the various gray pairs');
for i0 = 1:a
    for i = 1:b
        for j = 1:c
            fprintf('Error for gray pair [%g,%g,%g] is %g \n',...
                gray1(i0),gray2(i),gray3(j),levels((i0-1)*b*c+(i-1)*c+j,5));
        end
    end
end
