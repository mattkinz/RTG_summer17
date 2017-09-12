function [U,mu,beta,muf,betaf,muDbeta,sigma,...
    delta,gL,ind,out] = get_DIPSTV(p,q,r,Atb,scl,opts,k,b,wrap_shrink)

% initialize U, beta, mu
muf = opts.mu;       % final mu
betaf = opts.beta;     % final beta
beta = opts.beta0;
mu = opts.mu0;
muDbeta = mu/beta;





% initialize D^T sigma + A^T delta
gL = zeros(p*q*r,1);



% declare out variables
out.rel_chg_inn = [];  out.rel_chg_out = []; out.rel_lam2 = [];
out.f = [];        % record values of augmented Lagrangian fnc
out.cnt = [];      % record # of back tracking
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.tau = []; out.alpha = []; out.C = []; out.mu = [];
out.DU = [];


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

if round(k) == k & opts.levels == 1
    [D,Dt] = FD3D_dips(k,p,q,r);
elseif opts.levels == 1
    [D,Dt] = FFD3D(k,p,q,r); 
else
    [D,Dt] = FD3D_multiscale(k,p,q,r,opts.levels);
end

% if reweighted TV, override everything else
if opts.reweighted_TV
    [D,Dt] = FD3D_weighted(k,p,q,r,opts.coef);
end


% Check that Dt is the true adjoint of D
[flg,~,~] = check_D_Dt(D,Dt,[p,q,r]);
if ~flg
    error('Sparse domain transforms do not appear consistent');
end

% initialize multiplers
if isfield(opts,'sigma')
    sigma = opts.sigma;
else    
    sigma = D(zeros(p,q,r)); 
end
if isfield(opts,'delta')
    delta = opts.delta;
else
    delta = zeros(length(b),1);
end


if ~wrap_shrink
    ind = get_ind(k,p,q,r);
else
    ind=[];
end