function [U,mu,beta,beta2,muf,betaf,sigma,...
    delta,xi,gL,ind,out] = get_init_const(p,q,r,Atb,scl,opts,k,b,wrap_shrink)



% initialize U, beta, mu
muf = opts.mu;       % final mu
betaf = opts.beta;     % final beta
beta = opts.beta0;
mu = opts.mu0;
beta2 = beta;
%muDbeta = mu/beta;





% initialize D^T sigma + A^T delta
gL = zeros(p*q*r,1);



% declare out variables
out.rel_chg_inn = [];  out.rel_chg_out = []; out.rel_lam2 = [];
out.f = [];        % record values of augmented Lagrangian fnc
out.cnt = [];      % record # of back tracking
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.lam6 = []; out.lam7 = [];
out.tau = []; out.alpha = []; out.C = []; out.mu = [];
out.DU = [];

if opts.store_soln
    if r~=1
        warning('Solution is not being stored since r~=1')
        opts.store_soln = false;
    else
        out.Uall = zeros(p,q,opts.inner_iter*opts.outer_iter);
    end
end

% initialize U
[mm,nn,rr] = size(opts.init);
if max([mm,nn,rr]) == 1
    switch opts.init
        case 0, U = zeros(p,q,r);
        case 1, U = reshape(Atb,p,q,r);
    end
else
    if mm ~= p || nn ~= q || rr ~= r
        fprintf('Input initial guess has incompatible size! Switch to the default initial guess. \n');
        U = reshape(Atb,p,q,r);
    else
        U = opts.init*scl;
    end
end

global D Dt
if opts.smooth_phase
    if round(k) == k
        [D,Dt] = FD3D(k,p,q,r);
    else
        [D,Dt] = FFD3D(k,p,q,r); 
    end
else
    if sum(sum(sum(abs(opts.phase_angles))))==0
        opts.phase_angles = exp(-1i*reshape(Atb,p,q,r));
    else
        opts.phase_angles = exp(-1i*opts.phase_angles);
    end
    if round(k) == k
        [D,Dt] = FD3D_complex(k,p,q,r,opts.phase_angles);
    else
        [D,Dt] = FFD3D_complex(k,p,q,r,opts.phase_angles);
    end
end

if opts.adaptive
    [D,Dt] = FD3D_complex_reweight(k,p,q,r,opts.init);
end

% Check that Dt is the true adjoint of D
[flg,~,~] = check_D_Dt(D,Dt,[p,q,r]);
if ~flg
    error('Sparse domain transforms do not appear consistent');
end

% initialize multiplers
sigma = D(zeros(p,q,r));                       
delta = zeros(length(b),1);
xi = zeros(p*q*r,1);

if ~wrap_shrink
    ind = get_ind(k,p,q,r);
else
    ind=[];
end