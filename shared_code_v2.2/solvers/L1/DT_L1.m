function [U, out] = DT_L1(A,b,n,opts)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 05/17/2016


% This code has been modified to solve l1 penalty problems with the
% polynomial annihilation transform.  Several small bugs and notation
% changes have been made as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Problem Description       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original motivation to find:

%               min_f { mu/2 ||Af - b||_2^2 + ||D^k f||_1 }

% where D^k is kth order finite difference
% The problem is modified using variable splitting
% and this algorithm solves: 

%      min_{f,w} {mu/2 ||Af - b||_2^2 + beta/2 ||D^k f - w ||_2^2 
%               ||w||_1 - (delta , Af - b ) - (sigma , D^k f - w) }

% delta and sigma are Lagrange multipliers
% Algorithm uses alternating direction minimization over f and w.



% This algorithm was originally authored by Chengbo Li at Rice University.
% original code and description can be found here: 
% http://www.caam.rice.edu/~optimization/L1/TVAL3/

% Inputs: 
%   A: matrix operator as either a matrix or function handle
%   b: data values in vector form
%   p,q,r: signal dimensions
%   opts: structer containing opts, see function check_PA_opts.m for these


% Outputs:
%   U: reconstructed signal
%   out: output numerics

if numel(n)<3
    n(3) = 1;
end

p = n(1);q = n(2);r = n(3);


% get and check opts
opts = check_HOTV_opts(opts);


% mark important constants
tol_inn = opts.tol_inn;
tol_out = opts.tol_out;
k = opts.order;
n = p*q*r;
wrap_shrink = opts.wrap_shrink;
if round(k)~=k
    wrap_shrink = true;
end



% unify implementation of A
if ~isa(A,'function_handle')
    A = @(u,mode) f_handleA(A,u,mode);
end

% check scaling A
if opts.scale_A
    [A,b] = ScaleA(n,A,b);
end

% check scaling b
scl = 1;
if opts.scale_b
    [b,scl] = Scaleb(b);
end

% check for maximum constraint value
if opts.max_c
    max_v = opts.max_v*scl;
end


% calculate A'*b
Atb = A(b,2);


% initialize everything else
global D Dt
[U,mu,beta,muf,betaf,muDbeta,sigma,delta,gL,ind,out] ...
    = get_init(p,q,r,Atb,scl,opts,k,b,wrap_shrink);    % U: p*q

grays = sort(opts.grays)*scl;
ng = numel(grays);

nrmb = norm(b);
Upout = U;
Uc = D(U);
for i = 1:ng
    Uc(:,i+3) = Uc(:,i+3)-grays(i);
end


% first shrinkage step
W = max(abs(Uc) - 1/beta, 0).*sign(Uc);
% reset edge values if not using periodic regularization
if ~wrap_shrink, W(ind)=Uc(ind); end

lam1 = sum(sum(sum(abs(W))));

% gA and gD are the gradients of ||Au-b||^2 and ||Du-w||^2, respectively
% i.e. g = A'(Au-b), gD = D'(Du-w)
[lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
    lam1,beta,mu,A,b,Atb,sigma,delta);


% compute gradient
g = gD + muDbeta*gA - gL;


out.f = [out.f; f]; 
out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; 
out.lam3 = [out.lam3; lam3];out.lam4 = [out.lam4; lam4]; 
out.lam5 = [out.lam5; lam5];out.mu = [out.mu; mu];

        
for ii = 1:opts.outer_iter
    if opts.disp
            fprintf('    Beginning outer iteration #%d\n',ii);
            fprintf('    mu = %d , beta = %d , order = %g\n',mu,beta,k);
            fprintf('iter    ||w||_1    ||Du - w||^2  ||Au - b||^2\n');
    end
        
    %initialize the constants
    gam = opts.gam; Q = 1; fp = f;
    
    
    for jj = 1:opts.inner_iter
        % compute step length, tau
        if jj~=1
            % BB-like step length
            dgA = gA - gAp;   
            dgD = gD - gDp;                    
            ss = uup'*uup;                      
            sy = uup'*(dgD + muDbeta*dgA);       
            tau = abs(ss/max(sy,eps));          
        else
            % do Steepest Descent at the 1st ieration
            gc = D(reshape(g,p,q,r));       
            dDd = sum(sum(sum(gc.*conj(gc))));
            Ag = A(g,1);
            tau = abs((g'*g)/(dDd + muDbeta*(Ag')*Ag));
        end

        % keep previous values for backtracking & computing next tau
        Up = U; gAp = gA; gDp = gD; Aup = Au; 
        Ucp = Uc; %DtsAtdp =  DtsAtd;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ONE-STEP GRADIENT DESCENT %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        U = U(:) - tau*g;
        % projected gradient method for inequality constraints
        if opts.nonneg
            U = max(real(U),0);
        elseif opts.isreal
            U = real(U);
        end
        if opts.max_c
            U = min(U,max_v);
        end
        U = reshape(U,p,q,r);
        Uc = D(U);
        for i = 1:ng
            Uc(:,i+3) = Uc(:,i+3)-grays(i);
        end

        [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
            lam1,beta,mu,A,b,Atb,sigma,delta);

        % Nonmonotone Line Search Back tracking
        % Unew = Up + alpha*(U - Up)
        % f should be decreasing, if not, then the algorithm moves U
        % back in the direction of the previous solution
        alpha = 1;
        du = U - Up;
        const = 1e-5*beta*(g'*g*tau);
        cnt = 0; flg = true;
        %{
        while f > fp - alpha*const
            if cnt <5
                if flg
                    dgA = gA - gAp;
                    dgD = gD - gDp;
                    dAu = Au - Aup;
                    dUc = Uc - Ucp;
                    flg = false;
                end
                % shrink alpha
                alpha = alpha*opts.gamma;
                % U is moved back toward Up, in particular: 
                % U = alpha*U +(1-alpha)Up;
                % all other values are updated accordingly
                [U,lam2,lam3,lam4,lam5,f,Uc,Au,gA,gD] = back_up(p,q,r,...
                    lam1,alpha,beta,mu,Up,du,gAp,dgA,gDp,dgD,Aup,dAu,W,...
                    Ucp,dUc,b,sigma,delta);
                cnt = cnt + 1;
            else
                
                % shrink gam
                gam = opts.rate_gam*gam;

                % give up and take Steepest Descent step
                if (opts.disp > 0) && (mod(jj,opts.disp) == 0)
                    disp('    count of back tracking attains 5 ');
                end

                % compute step length, tau
                gc = D(reshape(g,p,q,r));
                dDd = sum(sum(sum(...
                    gc.*conj(gc))));
                Ag = A(g,1);
                tau = abs((g'*g)/(dDd + muDbeta*(Ag')*Ag));
                %update
                U = Up(:) - tau*g;
                % projected gradient method for inequality constraints
                if opts.nonneg
                    U = max(real(U),0);
                elseif opts.isreal
                    U = real(U);
                end
                
                U = reshape(U,p,q,r);
                Uc = D(U);
                Uc(:,end) = Uc(:,end) - scl;
                % shrinkage
                Ucbar = Uc - sigma/beta;
                W = max(abs(Ucbar) - 1/beta, 0).*sign(Ucbar);
                % reset edge values if not using periodic regularization
                if ~wrap_shrink, W(ind)=Uc(ind); end

                lam1 = sum(sum(sum(abs(W))));
                [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
                    lam1,beta,mu,A,b,Atb,sigma,delta);
                alpha = 0; % remark the failure of back tracking
                break;
            end
            
        end
        %}


        % if back tracking is successful, then recompute
        if alpha ~= 0
            Ucbar = Uc - sigma/beta;
            W = max(abs(Ucbar) - 1/beta, 0).*sign(Ucbar);
            % reset edge values if not using periodic regularization
            if ~wrap_shrink, W(ind)=Uc(ind); end
            % update parameters related to Wx, Wy
            [lam1,lam2,lam4,f,gD] = update_W(beta,...
                W,Uc,sigma,lam1,lam2,lam4,f);
        end

        % update reference value
        Qp = Q; Q = gam*Qp + 1; fp = (gam*Qp*fp + f)/Q;
        uup = U - Up; uup = uup(:);           % uup: pqr

        out.f = [out.f; f]; out.C = [out.C; fp]; out.cnt = [out.cnt;cnt];
        out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
        out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];
        out.tau = [out.tau; tau]; out.alpha = [out.alpha; alpha];out.mu = [out.mu; mu];

        if (opts.disp > 0) && (mod(ii,opts.disp) == 0)
            prnt_format = '%3.0f %10.5g %12.5g %13.5g\n';
            fprintf(prnt_format, jj,lam1,lam2,lam3);
        end


        % recompute gradient
        g = gD + muDbeta*gA - gL;
        
        % move to next outer iteration and update multipliers if relative
        % change is less than tolerance
        if (norm(g) < tol_inn), break; end;
        
        
    end
    % end of inner loop
    
    
    RelChgOut = norm(U(:)-Upout(:))/norm(Upout(:));
    out.reer = [out.reer; RelChgOut];
    Upout = U;

    % stop if already reached optimal solution
    if RelChgOut < tol_out || sqrt(lam3(end))/nrmb<opts.min_l2_error
        break;
    end

    % update multipliers
    [sigma,delta,lam4,lam5] = update_mlp(beta,mu, ...
        W,Uc,Au,b,sigma,delta);
    if ~opts.data_mlp, delta(:) = 0; lam5 = 0; end


    % update penality parameters for continuation scheme
    %beta0 = beta;
    beta = min(betaf, beta*opts.rate_ctn);
    mu = min(muf, mu*opts.rate_ctn);
    muDbeta = mu/beta;

    % update function value, gradient, and relavent constant
    f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
    %gL = -(beta0/beta)*g;     % DtsAtd should be divided by new beta  
    gL = 1/beta*(Dt(sigma) + A(delta,2));
    % gradient, divided by beta
    g = gD + muDbeta*gA - gL;

        
        
     
        

end

if sqrt(lam3(end))/nrmb<opts.min_l2_error
    fprintf('\nREACHED OPTIMAL L2 ERROR!!!\n\n');
end
out.total_iter = jj*ii;
out.final_error = norm(A(U(:),1)-b)/nrmb;
out.final_wl1 = lam1(end);
out.final_Du_w = lam2(end);
final_disp(out,opts);
            
% rescale U
U = U/scl;





function [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
    lam1,beta,mu,A,b,Atb,sigma,delta)
global Dt

% A*u
Au = A(U(:),1);

% gA 
gA = A(Au,2) - Atb;



% lam2, ||Du - w||^2
V = Uc - W;

%lam2 = sum(sum(sum(Vx.*Vx + Vy.*Vy + Vz.*Vz)));
lam2 = sum(sum(sum(V.*conj(V))));

% gD = D'(Du-w)
gD = Dt(V);



% lam3
Aub = Au-b;
lam3 = norm(Aub)^2;

%lam4
lam4 = sum(sum(sum(sigma.*V)));

%lam5
lam5 = delta'*Aub;

% f
f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;



function [U,lam2,lam3,lam4,lam5,f,Uc,Au,gA,gD] = back_up(p,q,r,lam1,...
    alpha,beta,mu,Up,du,gAp,dgA,gDp,dgD,Aup,dAu,W,Ucp,dUc,...
    b,sigma,delta)

gA = gAp + alpha*dgA;
gD = gDp + alpha*dgD;
U = Up + alpha*reshape(du,p,q,r);
Au = Aup + alpha*dAu;
Uc = Ucp + alpha*dUc;


V = Uc - W;


lam2 = sum(sum(sum(V.*conj(V))));
Aub = Au-b;
lam3 = norm(Aub)^2;
lam4 = sum(sum(sum(sigma.*V)));
lam5 = delta'*Aub;
f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;



function [lam1,lam2,lam4,f,gD] = update_W(beta,...
    W,Uc,sigma,lam1,lam2,lam4,f)
global Dt

% update parameters because Wx, Wy were updated
tmpf = f -lam1 - beta/2*lam2 + lam4;
lam1 = sum(sum(sum(abs(W))));
V = Uc - W;

gD = Dt(V);
lam2 = sum(sum(sum(V.*conj(V))));
lam4 = sum(sum(sum(sigma.*V )));
f = tmpf +lam1 + beta/2*lam2 - lam4;



function [sigma,delta,lam4,lam5] = update_mlp(beta,mu, ...
    W,Uc,Au,b,sigma,delta)


V = Uc - W;
sigma = sigma - beta*V;
Aub = Au-b;
delta = delta - mu*Aub;

%tmpf = f + lam4 + lam5;
lam4 = sum(sum(sum(sigma.*V)));
lam5 = delta'*Aub;
%f = tmpf - lam4 - lam5;







function [U,mu,beta,muf,betaf,muDbeta,sigma,...
    delta,gL,ind,out] = get_init(p,q,r,Atb,scl,opts,k,b,wrap_shrink)

% initialize U, beta, mu
muf = opts.mu;       % final mu
betaf = opts.beta;     % final beta
beta = opts.beta0;
mu = opts.mu0;
muDbeta = mu/beta;





% initialize D^T sigma + A^T delta
gL = zeros(p*q*r,1);



% declare out variables
out.reer = [];     % record relative errors of outer iterations
out.f = [];        % record values of augmented Lagrangian fnc
out.cnt = [];      % record # of back tracking
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.tau = []; out.alpha = []; out.C = []; out.mu = [];


% initialize U
[mm,nn,rr] = size(opts.init);

if mm ~= p || nn ~= q || rr ~= r
    error('Input initial guess has incompatible size!');    
else
    U = opts.init*scl;
end

if ~isfield(opts,'grays')
    error('gray values not specified');
end

global D Dt
[D,Dt] = DT_transform(k,p,q,r,opts.init,opts.grays);

% Check that Dt is the true adjoint of D
[flg,~,~] = check_D_Dt(D,Dt,[p,q,r]);
if ~flg
    error('Sparse domain transforms do not appear consistent');
end

% initialize multiplers
sigma = D(zeros(p,q,r));                       
delta = zeros(length(b),1);


if ~wrap_shrink
    ind = get_ind(k,p,q,r);
else
    ind=[];
end



function [] = final_disp(out,opts)

fprintf('Number of total iterations is %d. \n',out.total_iter);
fprintf('Final error: %5.3f\n',out.final_error);
fprintf('Final ||w||_1: %5.3g\n',out.final_wl1);
fprintf('Final ||Du-w||^2: %5.3g\n',out.final_Du_w);

if opts.disp
    figure(2);
    plot(out.lam1);  %plot lam1, ||W||_1
    hold on;
    plot(out.lam3.*out.mu);  %plot lam3, mu||Au -f||^2
    plot(abs(out.f));   %plot f, the objective function
    plot(out.mu)
    plot(2:opts.inner_iter:max(size(out.f)),out.f(2:opts.inner_iter:end),'kx','Linewidth',2);
    legend({'||W||_1','mu*||Au - f||_2^2',...
        'objective function','mu','multipliers updated'},...
        'Location','northeast');
    xlabel('iteration');
    title('Values in each iteration');
    hold off;

end