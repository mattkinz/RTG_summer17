function [D,Dt] = get_D_Dt(k,p,q,r,opts,Atb)


if opts.smooth_phase
    if round(k) == k & opts.levels == 1
        [D,Dt] = FD3D(k,p,q,r);  % standard FD
    elseif opts.levels == 1
        [D,Dt] = FFD3D(k,p,q,r); % fractional FD
    elseif round(k)==k
        if r~=1
            [D,Dt] = FFD3D_multiscale(k,opts.levels,p,q,r);  % multiscale FD
        else
            [D,Dt] = FFD2D_multiscale(k,opts.levels,p,q);  % multiscale 2D FD
        end
    else
        if r~=1
            [D,Dt] = FD3D_multiscale(k,opts.levels,p,q,r); % multiscale fractional FD
        else
            [D,Dt] = FD2D_multiscale(k,opts.levels,p,q); % multiscale fractional FD
        end
    end
else   % complex FD's using phase estimation
    % set up phase angles
    if sum(sum(sum(abs(opts.phase_angles))))==0
        opts.phase_angles = exp(-1i*reshape(Atb,p,q,r)); 
    else
        opts.phase_angles = exp(-1i*opts.phase_angles);
    end
    if round(k) == k && opts.levels == 1
        [D,Dt] = FD3D_complex(k,p,q,r,opts.phase_angles); % standard complex FD
    elseif opts.levels == 1
        [D,Dt] = FFD3D_complex(k,p,q,r,opts.phase_angles); % fractional complex FD
    else
        [D,Dt] = FD3D_multiscale_complex(k,opts.levels,p,q,r,opts.phase_angles);  % multiscale complex FD 
    end
end

% "joint" sparsity operators, i.e. different order in the z direction
if isfield(opts,'orderz')
    if opts.smooth_phase
        if k~=round(k) || opts.orderz ~= round(opts.orderz) 
            [D,Dt] = FFD3D_joint(k,opts.orderz,p,q,r);  % fractional FD
        else
            [D,Dt] = FD3D_joint(k,opts.orderz,p,q,r); % standard FD
        end
    else        
        % complex operators
        if k~=round(k) || opts.orderz ~= round(opts.orderz) 
            [D,Dt] = FFD3D_complex_joint(k,opts.orderz,p,q,r,opts.phase_angles);
        else
            [D,Dt] = FD3D_complex_joint(k,opts.orderz,p,q,r,opts.phase_angles);
        end
    end
end

% if reweighted TV, override everything else
if opts.reweighted_TV
    [D,Dt] = FD3D_weighted(k,p,q,r,opts.coef);
end