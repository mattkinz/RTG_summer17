function [Uchunk,Ucheck,set,fixed,update_type,type,convergence] = ...
    DIRT_thresh(Uchunk,opt,i,c_tol,Ucheck)

    convergence=0;
    if i<=opt.region_dirt_iter
        [Uchunk,set,fixed]=thresh_tol(Uchunk,opt.grays,c_tol);
        update_type=opt.region_update_type;
        type=[opt.region_update_type,', region=',num2str(c_tol)];
    else
        %{
        if opt.bad_pixels
            Uchunk = threshold(Uchunk,thresh,grays);
            [set,fixed]=bad_pixels1(W,bb,Uchunk,rand(1)*.75,opt.angles);
            opt.update_type = opt.region_update_type;
            type = [opt.update_type,', bad'];

        else
        %}
        %if i <50
        %[Uchunk,set,fixed]=threshbdryt(Uchunk,opt.thresh,opt.grays);
        %else
            [Uchunk,set,fixed]=threshbdry(Uchunk,opt.thresh,opt.grays);
        %end
        update_type = opt.bdry_update_type;
        type=[opt.bdry_update_type,', bdry'];
        %end
        %%check for convergence of the solution
        if opt.convergence_criteria
            if nnz(Ucheck-Uchunk)/numel(Uchunk)<opt.convergence_tol
                fprintf('CONVERGENCE MET!!! \n\n');
                convergence=1;
            else
                convergence=0;
            end
            Ucheck = Uchunk;
        end    
    end
    
%throw in some random pixels if needed, probably useless
        if opt.rand_pixel_probability
            for j = 1:size(Uchunk,3)
                ss = size(fixed{j},1);
                rrr = randperm(ss);
                rrr = rrr(1:round(ss*opt.rand_pixel_probability))';
                set{j}=[set{j};fixed{j}(rrr)];
                fixed{j}(rrr)='';
            end
        end
        
        
end