function [opt,stack] = check_dips_options(opt,stack)



%Written by: Toby Sanders
%Computational & Applied Mathematics Department, Univ. of South Carolina
%Dept. of Math and Stat Sciences, Arizona State University
%02/26/2016




    if ~isfield(opt,'thresh') || ~isfield(opt,'grays')
        error('User has not specified thresholds and/or gray values')
    end
    if max(size(opt.thresh))~= max(size(opt.grays))-1
        error('Number of thresholds and grays is inconsistent');
    end
    for i = 1:max(size(opt.thresh))
        if opt.thresh(i)<opt.grays(i) || opt.thresh(i)>opt.grays(i+1)
            error('thresholds are not between the gray values');
        end
    end
    
    % set maximum value constraints for the L1 solver
    opt.poly.max_v = opt.grays(end);
    opt.poly.min_v = opt.grays(1);
    
           
    if ~isfield(opt,'rx')
        opt.rx=3;
        opt.sigmax=1;
    end
    
    if ~isfield(opt,'rz')
        opt.rz=3;
        opt.sigmaz=1.5;
    end
    
    if ~isfield(opt,'chunksize')
        opt.chunksize=50;
        opt.overlap=5;
    end
    
    
    if ~isfield(opt,'disp')
        opt.disp=true;
    end
    if ~isfield(opt,'disppic')
        opt.disppic=false;
    end

    
    if ~isfield(opt,'angles')
        error('angles have not been specified');
    elseif min(opt.angles)<-90 || max(opt.angles)>90
        error('angles should be in the interval [-90,90]');
    end
    
    
    if stack==0
        if ~isfield(opt,'data_type')
            error('the type of file the tilt series is saved is not specified');
        elseif strcmp(opt.data_type,'mat')
                stack=load(opt.data_name);
                stack = struct2cell(stack);
                if max(size(stack))>1
                    error('matlab file for the stack has multiple variables');
                end
                stack = stack{1};
        elseif strcmp(opt.data_type,'mrc')
            stack=ReadMRC(opt.mrcfilename);
            stack=double(stack);
        else
            error('opt.data_type should be "mat" or "mrc"');
        end
    end
    
    if ~isfield(opt,'resolution')
        opt.resolution=size(stack,1);
        fprintf('Reconstruction size not specified, setting to %i\n',size(stack,1));
    end
        
    if isfield(opt,'testslice')
        opt.startslice=opt.testslice;
        opt.endslice=opt.testslice;
    else
        if ~isfield(opt,'startslice')
            opt.startslice=1;
            fprintf('no startslice specified,\n reconstruction starting from beginning\n\n');
        elseif opt.startslice>size(stack,2) 
             opt.startslice=1;
             fprintf('input startslice was larger than slices in the stack,\n setting startslice to the beginning\n\n');
        end
         if ~isfield(opt,'endslice')
            opt.endslice=size(stack,2);
            fprintf('no endslice specified,\n reconstruction finishing at end of stack\n\n');
        elseif opt.endslice>size(stack,2)
            opt.endslice = size(stack,2);
            fprintf('input endslice was larger than slices in the stack,\n setting endslice to the end\n\n');
        end
    end
    
    if ~isfield(opt,'initialsolution')
        opt.initialsolution=false;
    end
    
    if ~isfield(opt,'model')
        opt.model = 'line';
    elseif ~ischar(opt.model)
        opt.model = 'line';
    elseif ~sum(strcmp(opt.model,{'line','strip'}))
        opt.model = 'line';
    end
    
    if ~isfield(opt,'bdry_dirt_iter')
        opt.bdry_dirt_iter=0;
    end
    if ~isfield(opt,'region_dirt_iter')
        opt.region_dirt_iter=0;
    end
    
    
    
    if ~isfield(opt,'t_tol')
        mm=opt.grays(end)-opt.grays(1);
        for i = 1:max(size(opt.grays))-1
            if opt.grays(i+1)-opt.grays(i)<mm
                mm = opt.grays(i+1)-opt.grays(i);
            end
        end
        opt.t_tol=mm/5;
        fprintf('region tolerance not specified \n, set to %d',mm/5);
    else
        mm=opt.grays(end)-opt.grays(1);
        for i = 1:max(size(opt.grays))-1
            if opt.grays(i+1)-opt.grays(i)<mm
                mm = opt.grays(i+1)-opt.grays(i);
            end
        end
        if mm/2<=opt.t_tol
            error('opt.t_tol is set too high');
        end
        if mm*1/4<=opt.t_tol
            fprintf('\n warning: opt.t_tol is very high\n');
        end
    end
    if ~isfield(opt,'t_delta')
        opt.t_delta = 1/4*opt.t_tol;
    end
    
    if ~isfield(opt,'t_epsilon')
        opt.t_epsilon=.1;
    end
    
    
    
    %PA solver options
    if ~isfield(opt.poly,'order')
        opt.poly.order=1;
        fprintf(['order of the finite difference operator not set\n',...
            'using default order 1\n']);
    end

    if ~isfield(opt.poly,'nonneg')
        opt.poly.nonneg = false;
    end
    
    if ~isfield(opt.poly,'mu')
        opt.poly.mu = 60*2^(opt.poly.order-1);
        fprintf('mu not set, using default\n');
    end
    
    
        
    
    

