function [opt,stack] = checkdirtoptions(opt,stack)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014

%Last updated on 4/30/2015.

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
    
    if ~isfield(opt,'bdry_dirt_iter')
        fprintf('\nboundary iterations not specified,\n using convergence criteria\n');
        opt.convergence_criteria=true;
    end
    
    if ~isfield(opt,'region_dirt_iter')
        opt.region_dirt_iter = 0;
    elseif opt.region_dirt_iter>0
        if ~isfield(opt,'t_tol')
            mm=opt.grays(end)-opt.grays(1);
            for i = 1:max(size(opt.grays))-1
                if opt.grays(i+1)-opt.grays(i)<mm
                    mm = opt.grays(i+1)-opt.grays(i);
                end
            end
            opt.t_tol=mm/4;
            fprintf('region tolerance not specified \n, set to %d',mm/4);
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
            opt.t_epsilon=.01;
        end
    
    end
    

    
    
    if ~isfield(opt,'convergence_criteria')
        opt.convergence_criteria=false;
        opt.convergence_tol=10^(-6);
        opt.total_iter = opt.bdry_dirt_iter + opt.region_dirt_iter;
    elseif ~opt.convergence_criteria
        opt.total_iter = opt.bdry_dirt_iter + opt.region_dirt_iter;
    else
        if ~isfield(opt,'total_iter')
            opt.total_iter=50;
        end
        if ~isfield(opt,'convergence_tol')
            opt.convergence_tol=10^(-6);
        end
    end
        

    
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
        opt.overlap=10;
    end
    
    
    
    
    
    if ~isfield(opt,'disp')
        opt.disp=true;
    end
    if ~isfield(opt,'disppic')
        opt.disppic=true;
    end
    
    if ~isfield(opt,'angles')
        error('angles have not been specified');
    elseif min(opt.angles)<-90 || max(opt.angles)>90
        error('angles should be in the interval (-90,90)');
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
            stack=ReadMRC(opt.data_name);
            stack=double(stack);
        else
            error('opt.filetype should be "mat" or "mrc"');
        end
    end
    
    
    
    if ~isfield(opt,'resolution')
        opt.resolution=size(stack,1);
        fprintf('Resolution not specified, setting to %i\n',size(stack,1));
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
    
    
    if ~isfield(opt,'region_update_type')
        opt.region_update_type = 'SIRT';
    elseif ~sum(strcmp(opt.region_update_type,{'CGLS','SIRT'}))
        opt.region_update_type = 'SIRT';
        fprintf('\n Region update type not recognized, set to SIRT\n');
    end
    
    if ~isfield(opt,'bdry_update_type')
        opt.bdry_update_type = 'CGLS';
    elseif ~sum(strcmp(opt.bdry_update_type,{'CGLS','SIRT'}))
        opt.bdry_update_type = 'CGLS';
        fprintf('\n Region update type not recognized, set to CGLS\n');
    end
    
    if ~isfield(opt,'inner_iter');
        opt.inner_iter = 20;
    elseif opt.inner_iter<10
        fprintf('\n warning: recommended that inner iterations at least 10\n');
    end
    
    if ~isfield(opt,'rand_pixel_probability')
        opt.rand_pixel_probability = 0;
    end
  
    if ~isfield(opt,'minc')
        opt.minc=true;
    end
    
    
end
    
    

    
        
    
    

