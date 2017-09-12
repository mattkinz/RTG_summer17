function [opt,stack] = checkpdirtoptions(opt,stack)



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014



    if ~isfield(opt,'thresh') || ~isfield(opt,'grays')
        error('User has not specified thresholds and/or gray values')
    end
    if max(size(opt.thresh))~= max(size(opt.grays))
        error('Number of thresholds and grays is inconsistent for PDIRT');
    end

    
    if ~isfield(opt,'bdry_dirt_iter')
        fprintf('\nboundary iterations not specified,\n using convergence criteria\n');
        opt.convergence_criteria=true;
    end
    
    if ~isfield(opt,'total_iter')
        opt.total_iter = 20;
        fprintf('\n opt.total_iter not specified, set to 20 iterations\n');
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
    if ~isfield(opt,'initial_solution')
        opt.initial_solution=false;
    end
    if ~isfield(opt,'model')
        opt.model = 'line';
    elseif ~ischar(opt.model)
        opt.model = 'line';
    elseif ~sum(strcmp(opt.model,{'line','strip'}))
        opt.model = 'line';
    end
    
    
    if ~isfield(opt,'update_type')
        opt.update_type = 'SIRT';
    elseif ~strcmp(opt.update_type,{'CGLS','SIRT'})
        opt.update_type = 'CGLS';
        fprintf('\n Update type not recognized, set to SIRT\n');
    elseif strcmp(opt.update_type,'CGLS')
        fprintf('\nWARNING!!! Recommended SIRT update for PDART\n');
    end
    
    if ~isfield(opt,'inner_iter');
        opt.inner_iter = 20;
    elseif opt.inner_iter<10
        fprintf('\n warning: recommended that inner iterations at least 10\n');
    end
  
    
end
    
    

    
        
    
    

