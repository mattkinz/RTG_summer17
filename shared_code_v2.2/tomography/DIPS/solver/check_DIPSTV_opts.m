function opts = check_DIPSTV_opts(opts)

% Written by Toby Sanders @ASU
% Department of Statistical & Mathematical Sciences
% 05/17/2016

% order of the finite difference operator
% Set to 1 for TV regularization
% Set to 0 for signal sparsity
% Set to >= 2 for higher order PA methods
% Noninteger values are also accepted for fractional finite difference
% scheme


wflg = zeros(10,1);
if ~isfield(opts,'order')
    fprintf('Order of finite difference set to 1 (TV regularization)\n\n');
    opts.order=1;
elseif opts.order<0
    error('opts.order should be at least 0');
end



% scale mu according to the order of PA transform.  Highly recommended.
if isfield(opts,'scale_mu')
    if ~islogical(opts.scale_mu);
        error('opts.scale_mu should be true or false.');
    end
else
    opts.scale_mu = false;
end



% mu is generally the most important parameter
% mu is mainly decided by noise level. Set mu big when b is noise-free
% whereas set mu small when b is very noisy.
if isfield(opts,'mu')
    if ~isscalar(opts.mu) || opts.mu <0
        error('opts.mu must be positive.');
    elseif opts.scale_mu
        if opts.mu < 24  || opts.mu > 201
            wflg(1) = 1;           
        end
    else
        if opts.mu*2^(1-opts.order) < 24 || opts.mu*2^(1-opts.order)> 201
            wflg(1) = 1;
        end
    end
else
    %default mu
    if opts.scale_mu
        opts.mu = 90;
    else
        opts.mu = 90*2^(opts.order-1);
    end
end



% initial mu for continuation scheme
% method may yield better convergence starting with smaller mu
if isfield(opts,'mu0')
    if ~isscalar(opts.mu0) || opts.mu0 <= 0
        error('opts.mu0 is should be a positive number which is no bigger than beta.');
    end
else
    opts.mu0 = opts.mu/4;  
end



% beta can also be scaled according to the order
% This is recommended true and improves the results slightly
if isfield(opts,'scale_beta')
    if ~islogical(opts.scale_beta);
        error('opts.scale_beta should be true or false.');
    end
else
    opts.scale_beta = false;
end


% coefficient for sparsifying operator
% setting beta = 2^5 usually works well
if isfield(opts,'beta')
    if ~isscalar(opts.beta) || opts.beta <0
        error('opts.beta must be positive.');
    else if opts.beta > 2^13 || opts.beta < 2^4
            wflg(2) = 1;            
        end
    end
else
    % default beta
    opts.beta = 2^5;
end



% initial beta
if isfield(opts,'beta0')
    if ~isscalar(opts.beta0) || opts.beta0 <= 0
        error('opts.beta0 is should be a positive number which is no bigger than beta.');
    end
else
    opts.beta0 = opts.beta; 
end



% option for periodic regularization
% Typically we do not want to apply the shrinkage at the boundaries
% In this case set wrap_shrink to false
if ~isfield(opts,'wrap_shrink')
    opts.wrap_shrink = false;
elseif opts.wrap_shrink
    % Warn the user!
    wflg(3) = 1;
end

% number of inner loop iterations
% inner_iter gives the number of iterations (alternating gradient decent 
% with shrinkage) for each set of Lagrangian multipliers
if isfield(opts,'inner_iter')
    if ~isscalar(opts.inner_iter) || opts.inner_iter <= 0
        error('opts.inner_iter should be a positive integer.');
    end
else
    opts.inner_iter = 10;
end

% Maximum number of outer iterations (Updates on Lagrangian multipliers)
if isfield(opts,'outer_iter')
    if ~isscalar(opts.outer_iter) || opts.outer_iter <= 0
        error('opts.outer_iter should be a positive integer.');
    end
else
    opts.outer_iter = 10;
end


% continuation parameter for mu and beta
% after each outer iteration, mu = min(muf , mu*opts.rate_ctn)
if isfield(opts,'rate_ctn')
    if ~isscalar(opts.rate_ctn) || opts.rate_ctn <= 1
        error('opts.rate_ctn is either not a scalar or no bigger than one.');
    end
else
    opts.rate_ctn = 1.5;
end





% min/max value constraint
% for discrete tomography the default is true
if ~isfield(opts,'max_c')
    opts.max_c = true;
end

if ~isfield(opts,'min_c')
    opts.min_c = true;
end


% Convergence based on the relative l2 error
if isfield(opts,'min_l2_error')
    if opts.min_l2_error<0 || opts.min_l2_error>=1
        error('opts.min_l2_error should be between 0 and 1');
    elseif opts.min_l2_error>.2
        wflg(4) = 1;
    end
else
    opts.min_l2_error=0;
end


% outer loop convergence tolerance
if isfield(opts,'tol_out')
    if ~isscalar(opts.tol_out) || opts.tol_out <= 0
        error('opts.tol_out should be a positive small number.');
    end
else
    opts.tol_out = 1e-3;
end;


% inner loop convergence tolerance
if isfield(opts,'tol_inn')
    if ~isscalar(opts.tol_inn) || opts.tol_inn <= 0
        error('opts.tol_inn should be a positive small number.');
    end
else
    opts.tol_inn = 3e-3;
end;










% if the user has an initial guess, store it in this option
if isfield(opts,'init')
    if numel(opts.init) ~= 1
        %fprintf('User has supplied opts.init as initial guess solution!!!\n');
    end
else
    opts.init = 1;
end





% Display option
if ~isfield(opts,'disp')
    opts.disp = false;
end



% The user has the option to use Lagrangian multipliers for the data term.
% Default is to use the Lagrangian multipliers.  There is no option for
% this for the sparsity term.  The sparsity multipliers should be used to
% enforce the constrained problem, Du = w.
if ~isfield(opts,'data_mlp')
    opts.data_mlp = opts.outer_iter;
elseif ~opts.data_mlp
    wflg(7) = 1;
end
    



% number of scalings in the finite difference operator
if ~isfield(opts,'levels')
    opts.levels = 1;
elseif round(opts.levels)~=opts.levels || opts.levels < 1
    error('opts.levels should be a positive integer');
end

% if adaptive is true, then the method is using reweighted FD transform
% the new weights should be put into opts.coef, which is a 3x1 cell
if ~isfield(opts,'reweighted_TV')
    opts.reweighted_TV=false;
elseif opts.reweighted_TV == true & ~isfield(opts,'coef')
    error('Reweighted norm was set true without specifying the weights');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The remaining parameters are for the gradient decent and backtracking %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults for these parameters is recommended













% gamma is for backtracking
if isfield(opts,'gamma')
    if ~isscalar(opts.gamma) || opts.gamma <= 0 || opts.gamma > 1
        error('opts.gamma should be a scalar between 0 and 1.');
    end
else
    opts.gamma = .6;
end



% Control the degree of nonmonotonicity. 0 corresponds to monotone line search.
% The best convergence is obtained by using values closer to 1 when the iterates
% are far from the optimum, and using values closer to 0 when near an optimum.
if isfield(opts,'gam')
    if ~isscalar(opts.gam) || opts.gam <= 0 || opts.gam > 1
        error('opts.gam should be a scalar between 0 and 1.');
    end
else
    opts.gam = .9995;
end



% shrinkage rate of gam
if isfield(opts,'rate_gam')
    if ~isscalar(opts.rate_gam) || opts.rate_gam <= 0 || opts.rate_gam > 1
        error('opts.rate_gam should be a scalar between 0 and 1.');
    end
else
    opts.rate_gam = .9;
end



% tau is the step length in the gradient decent
if isfield(opts,'tau')
    if ~isscalar(opts.tau) || opts.tau <= 0
        error('opts.tau is not positive scalar.');
    end
else
    opts.tau = 1.8;
end
opts.StpCr=0;



% scale mu and beta here
ko = opts.order;
if opts.scale_mu && ko ~= 0
    opts.mu = opts.mu*2^(ko-1)*opts.levels;
    opts.mu0 = opts.mu0*2^(ko-1)*opts.levels;    
end

if opts.scale_beta && ko ~= 0
   % opts.beta = opts.beta*...
   %     factorial(ko-1)^2.45/factorial(2*(ko-1));
   % opts.beta0 = opts.beta0*...
   %     factorial(ko-1)^2.45/factorial(2*(ko-1));
    opts.beta = opts.beta*2^(1-ko);%*opts.levels;
    opts.beta0 = opts.beta0*2^(1-ko);%*opts.levels;
end

%{
if wflg1, fprintf('opts.mu is not within optimal range\n'); end;
if wflg2, fprintf('opts.beta is not within optimal range\n'); end;
if wflg3, fprintf('Using periodic regularization\n'); end;
if wflg4, fprintf('opts.min_l2_error is large'); end;
if wflg5, fprintf('scale_A set to false, this is not recommended'); end;
if wflg6, fprintf('scale_b set to false, this is not recommended'); end;
if wflg(7), warning('Lagrangian multiplier is not being used for data
enforcement'); end;
if wflg(8),         warning('signal is assumed to be complex with a smooth
phase'); end;
wflg(9), warning('signal set to be real but smooth_phase set false.');

wflg(10), warning('correction phase angles not specified, using Atb');
%}

msgs = ...
    {'-opts.mu is not within optimal range';
    '-opts.beta is not within optimal range';
    '-Using periodic regularization';
    '-opts.min_l2_error is large';
    '-scale_A set to false (not recommended)';
    '-scale_b set to false (not recommended)';
    '-Lagrangian multiplier is not being used to encourage the constrained problem';
    '-signal is assumed to be complex with a smooth phase';
    '-signal set to be real but smooth_phase set false'
    '-estimate phase angles not specified, using Atb'};



    
if sum(wflg)
    fprintf('WARNINGS:\n');
    disp(msgs(find(wflg)));
    %msgbox(msgs(find(wflg)),'PA3D: WARNINGS');
end
    








