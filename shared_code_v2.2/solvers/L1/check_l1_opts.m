function opts = check_l1_opts(opts)


% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


if ~isfield(opts,'D')
    error('Sparse transform operators not defined');
elseif ~isfield(opts,'Dt')
    error('Transpose sparse transform not defined');
end

wflg = zeros(10,1);
if ~isfield(opts,'order')
    fprintf('Order of finite difference set to 1 (TV regularization)\n\n');
    opts.order=1;
elseif opts.order<0
    error('opts.order should be at least 0');
end


% mu is generally the most important parameter
% mu is mainly decided by noise level. Set mu big when b is noise-free
% whereas set mu small when b is very noisy.
if isfield(opts,'mu')
    if ~isscalar(opts.mu) || opts.mu <0
        error('opts.mu must be positive.');
    end
else
    %default mu
    opts.mu = 90;
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
opts.wrap_shrink = true;

% number of inner loop iterations
% inner_iter gives the number of iterations (alternating gradient decent 
% with shrinkage) for each set of Lagrangian multipliers
if isfield(opts,'inner_iter')
    if ~isscalar(opts.inner_iter) || opts.inner_iter <= 0
        error('opts.inner_iter should be a positive integer.');
    end
else
    opts.inner_iter = 15;
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





% Nonnegativity constraint
if isfield(opts,'nonneg')
    if ~islogical(opts.nonneg)
        error('opts.nonneg should be true or false.');
    end
else
    opts.nonneg = false;
end

% Maximum value constraint
if ~isfield(opts,'max_c')
    opts.max_c = false;
elseif opts.max_c
    if ~isfield(opts,'max_v')
        error('maximum constraint (max_c) was set to true without specifying the maximum value (max_v)');
    end
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
    opts.tol_out = 1.e-3;
end;





% inner loop convergence tolerance
if isfield(opts,'tol_inn')
    if ~isscalar(opts.tol_inn) || opts.tol_inn <= 0
        error('opts.tol_inn should be a positive small number.');
    end
else
    opts.tol_inn = 1.e-3;
end;










% if the user has an initial guess, store it in this option
if isfield(opts,'init')
    if numel(opts.init) ~= 1
        fprintf('User has supplied opts.init as initial guess solution!!!\n');
    end
else
    opts.init = 1;
end





% Display option
if ~isfield(opts,'disp')
    opts.disp = false;
end

% Display figures options
if ~isfield(opts,'disp_fig')
    opts.disp_fig = opts.disp;
end




% scaling of the operator A and vector b. Scaling HIGHLY recommended so that
% consistent values for mu and beta may be used independent of the problem
% A is scaled so that ||A||_2 = 1.
if isfield(opts,'scale_A')
    if ~islogical(opts.scale_A)
        error('opts.scale_A should be true or false.');
    end
    if ~opts.scale_A
        wflg(5) = 1;        
    end
else
    opts.scale_A = true;
end


if isfield(opts,'scale_b')
    if ~islogical(opts.scale_b)
        error('opts.scale_b should be true or false.');
    end
    if ~opts.scale_b
        wflg(6) = 1;
    end
    
else
    opts.scale_b = true;
end




% option to store the solution at each iterate
% generally this should be false since it requires significant memory
if isfield(opts,'store_soln')
    if ~islogical(opts.store_soln)
        error('opts.store_soln should be true or false.');
    end
else
    opts.store_soln = false;
end




% The user has the option to use Lagrangian multipliers for the data term.
% Default is to use the Lagrangian multipliers.  There is no option for
% this for the sparsity term.  The sparsity multipliers should be used to
% enforce the constrained problem, Du = w.
if ~isfield(opts,'data_mlp')
    opts.data_mlp = opts.outer_iter;
elseif ~opts.data_mlp
    wflg(7) = 1;
elseif islogical(opts.data_mlp)
    if opts.data_mlp
        opts.data_mlp = opts.outer_iter;
    end
end
    


%if the data is complex, such as fourier data, but the signal is real, set
%this option to true to search for a real solution
if isfield(opts,'isreal')
    if ~islogical(opts.isreal)
        error('opts.isreal should be true or false.');
    end
    if ~opts.isreal && opts.smooth_phase
        wflg(8) = 1;
    end
else
    opts.isreal = false;
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




% if the signal is not necessarily real, its absolute value may be smooth
% but the phase at each pixel may be totally random.  In this case, set 
% smooth_phase to false, and specify the estimated phase angle (in radians) 
% at each pixel into opt.phase_angles
if ~isfield(opts,'smooth_phase')
    opts.smooth_phase = true;
elseif ~opts.smooth_phase
    if opts.isreal
        wflg(9) = 1;       
        %opts.smooth_phase = true;
    end
end


if ~isfield(opts,'phase_angles')
    opts.phase_angles = false;
end

if ~opts.smooth_phase && sum(sum(sum(opts.phase_angles)))==0
    wflg(10) = 1;    
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



    



if ~opts.isreal & ~opts.smooth_phase
    fprintf('***********************************\n')
    fprintf('*    PHASE ESTIMATED FD scheme    *\n');
    fprintf('***********************************\n');
end

if sum(wflg)
    fprintf('WARNINGS:\n');
    disp(msgs(find(wflg)));
    %msgbox(msgs(find(wflg)),'PA3D: WARNINGS');
end
    


