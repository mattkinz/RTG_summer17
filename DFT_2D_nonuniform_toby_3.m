%This script attempts to take the squarepulse function and make it more
%complicated by defining other types of functions across the square unit
%pulse

clear
clc
  

N =100;
%need to keep q = r = 1
q = 1;  %frequencies
r = 1;  %space

x = -pi + (2*pi/(2*r*N+1))*(0:2*r*N);  %physical space is a square
y = x;

%Define all of the parameters for each rectangle of the square
%A*squarepulse, B1*x, B2*y, C1*x^2, C2*y^2
%background_left
Abl = 16 - .25; B1bl = 8; B2bl = 0; C1bl = 1; C2bl = 0;
%background_right
Abr = 16 - .25; B1br = -8; B2br = 0; C1br = 1; C2br = 0;
%background_top
Abt = 2*pi - .75; B1bt = 0; B2bt = -1; C1bt = 0; C2bt = 0;
%background_bottom
Abb = -.75; B1bb = 0; B2bb = -1; C1bb = 0; C2bb = 0;
%upper_box
Aub = 4.5; B1ub = 0; B2ub = -6; C1ub = .5; C2ub = 2;
%lower_box
Alb = 4.5; B1lb = 0; B2lb = 6; C1lb = -.5; C2lb = 2;
%left_center_box
Alc = 1.5^2 + 1; B1lc = 3; B2lc = 0; C1lc = 1; C2lc = 0;
%middle_center_box
Amc = 0; B1mc = 0; B2mc = 0; C1mc = 8; C2mc = 0;
%right_center_box
Arc = .5; B1rc = 3; B2rc = 0; C1rc = 0; C2rc = 0;


%% Function values, calls functions that are stored in directory F and uses them
% to build the function F defined over a square`
F = zeros(numel(x));
for k = 1:numel(x)
	for j = 1:numel(y)
		if x(k) <= -2.5
            F(j,k) = f_background_left_piecewise(x(k),y(j),Abl,B1bl,B2bl,C1bl,C2bl);
		elseif x(k) > -2.5 && x(k) <= 2.5
			if y(j) > 2.5
				F(j,k) = f_background_top_piecewise(x(k),y(j),Abt,B1bt,B2bt,C1bt,C2bt);
			elseif y(j) <= 2.5 && y(j) > .5
				F(j,k) = f_upper_box_piecewise(x(k),y(j),Aub,B1ub,B2ub,C1ub,C2ub);
			elseif y(j) <= .5 && y(j) > -.5
				if x(k) <= -.5
					F(j,k) = f_left_center_box_piecewise(x(k),y(j),Alc,B1lc,B2lc,C1lc,C2lc);
				elseif x(k) <= .5
					F(j,k) = f_middle_center_box_piecewise(x(k),y(j),Amc,B1mc,B2mc,C1mc,C2mc);
				else
					F(j,k) = f_right_center_box_piecewise(x(k),y(j),Arc,B1rc,B2rc,C1rc,C2rc);
                end
			elseif y(j) <= -.5 && y(j) > -2.5
				F(j,k) = f_lower_box_piecewise(x(k),y(j),Alb,B1lb,B2lb,C1lb,C2lb);
			else
				F(j,k) = f_background_bottom_piecewise(x(k),y(j),Abb,B1bb,B2bb,C1bb,C2bb);
            end
        else
			F(j,k) = f_background_right_piecewise(x(k),y(j),Abr,B1br,B2br,C1br,C2br);
		end
	end
end
figure(1)
mesh(x,y,F);colorbar;
title('F(x,y)','interpreter','latex','fontsize',16)
xlabel('x','interpreter','latex','fontsize',16)
ylabel('y','interpreter','latex','fontsize',16)
		
		

%----------------------------------------------------------------------
%This section generates the pattern for which we will sovle for the fourier coefficients 
%using the exact fourier transform (integral) for our simple functions 
%------------------------------------------------------------------------

%% Frequency values, gridded
%the Fhat functions take vectors of frequency samples as input

% %integer spaced frequencies
%  wx = (-q*N:q*N)';
%  wy = (-q*N:q*N)';
%  [WX,WY] = meshgrid(wx,wy);
%  kx = WX(:);
%  ky = WY(:);
%-------------------------------------------------------------------------
%randomly generated frequencies about each integer point
% epsilon = .15;
% wx = rand_on_integers(epsilon,q,N).';
% wy = rand_on_integers(epsilon,q,N).';   %these have been turned into column vectors
% [WX,WY] = meshgrid(wx,wy);
% kx = WX(:);
% ky = WY(:);
%--------------------------------------------------------------------------
%frequencies distrubuted by the inverse square
% wx = (power_frequencies(q,N,2)).';
% wy = (power_frequencies(q,N,2)).';
% [WX,WY] = meshgrid(wx,wy);
% kx = WX(:);
% ky = WY(:);
%-------------------------------------------------------------------------
% %random points in k-space
% K = N*(2*rand((2*N+1)^2,2)-1);
% kx = K(:,1);
% ky = K(:,2);
%-------------------------------------------------------------------------
% %Spiral sampling
z = .05;
w = ((N+1)^2) ;
theta = linspace(0.01,N/z,w).';
kx = z*theta.*cos(theta);
ky = z*theta.*sin(theta);
figure(3)
xlim([-N,N])
ylim([-N,N])
plot(kx,ky,'b.')
%-------------------------------------------------------------------------
% Rosette
% taq = 0.05;
% theta = linspace(0,0.05,round(taq/(5e-6)));
% o1 = 3576.7;
% o2 = 3967.9;
% kx = (N*cos(o2*theta).*cos(o1*theta)).';
% ky = (N*cos(o2*theta).*sin(o1*theta)).';
% figure(1)
% xlim([-N,N])
% ylim([-N,N])
% plot(kx,ky,'b.')
%-------------------------------------------------------------------------
% %Gaussian points in k-space
% K = N/3*randn((2*N+1)^2,2);
% kx = K(:,1);
% ky = K(:,2);
%% Fhat vector 
%Fhat = zeros(sqrt(numel(kx)));

%============================================================================
%This section is where we generate the fourier coefficients. We calculate the 
%coefficients by solving the fourier integral over each section. This is done
%in the functions defined in the directory Fhat
%===========================================================================
Fhat = fhat_background_left_piecewise(kx,ky,Abl,B1bl,B2bl,C1bl,C2bl);
Fhat = Fhat + fhat_background_right_piecewise(kx,ky,Abr,B1br,B2br,C1br,C2br);
Fhat = Fhat + fhat_background_top_piecewise(kx,ky,Abt,B1bt,B2bt,C1bt,C2bt);
Fhat = Fhat + fhat_background_bottom_piecewise(kx,ky,Abb,B1bb,B2bb,C1bb,C2bb);
Fhat = Fhat + fhat_upper_box_piecewise(kx,ky,Aub,B1ub,B2ub,C1ub,C2ub);
Fhat = Fhat + fhat_lower_box_piecewise(kx,ky,Alb,B1lb,B2lb,C1lb,C2lb);
Fhat = Fhat + fhat_left_center_box_piecewise(kx,ky,Alc,B1lc,B2lc,C1lc,C2lc);
Fhat = Fhat + fhat_middle_center_box_piecewise(kx,ky,Amc,B1mc,B2mc,C1mc,C2mc);
Fhat = Fhat + fhat_right_center_box_piecewise(kx,ky,Arc,B1rc,B2rc,C1rc,C2rc);

%add random noise to my fourier coefficients

sigma= 1e-3*norm(Fhat(:),inf);
for j = 1:sqrt(numel(kx))
    Fhat(:,j) = Fhat(:,j) + sigma*(randn(sqrt(numel(kx)),1)) + sigma*(randn(sqrt(numel(kx)),1))*1i; 
end

%% Build the Fourier matrix, compuationally expensive and impractical to invert 

%E = exp(-1i*kx(:)*xx(:)'-1i*ky(:)*yy(:)')/(2*N+1)^2;


%% NUFFT utilities: utilizing the NUFFT programs in the directory irt
%build structure for nufft and nufft_adj
J = 5;		% interpolation neighborhood
K1 = 2*(2*N+1);	% two-times oversampling
%
st = nufft_init([kx,ky]*2*pi/(2*N+1), [2*N+1,2*N+1], [J J], [K1 K1],[(2*N+1)/2,(2*N+1)/2]);
%
%% Toby Sanders code, using the regularization programs in the directory shared_code_v2.2

 % options for l1 optimization algorithm
 clear pat;clear opt;
 pat.nonneg = false;
 pat.levels = 3;
 pat.order = 2;
 %pat.L1type = 'isotropic';
 pat.inner_iter = 10;
 pat.outer_iter = 17;
 pat.tol_inn = 1e-9;
 pat.tol_out = 1e-9;
 pat.disp = true;
 pat.isreal = true;
 pat.data_mlp = false;
 
mu = 2500;
%D = gallery('tridiag',zeros(numel(x)-1,1),-1*ones(numel(x),1),ones(numel(x)-1,1));


pat.mu = mu;
rec1 = (2*N+1)^2 * HOTV3D(@(u,mode) f_handleA_matt_2(u,mode,st,N),Fhat(:),[2*N+1,2*N+1,1],pat);
norm_2 = norm(F - rec1,'fro');
norm_1 = norm(F - rec1,1);
%norm_1_reg = norm(D*rec1,1);

figure(2)
mesh(x,y,rec1);colorbar;
title('F reconstructed, spiral, order 2, levels 3, mu = 2500','interpreter','latex','fontsize',16)
xlabel('x','interpreter','latex','fontsize',16)
ylabel('y','interpreter','latex','fontsize',16)

 
 % compute reconstructions for each order and each method (HOTV,
 % MHOTV and Daub wavelets)
%  for k = 1     % loop over the levels
%      
%      % HOTV reconstruction
%      pat.levels = k;
%      pat.mu = mu;
%      rec1(:,:,k) = (2*N+1)^2 * HOTV3D(@(u,mode) f_handleA_matt_2(u,mode,st,N),Fhat(:),[2*N+1,2*N+1,1],pat);
% end
 
%  
% figure(2)
% set(gcf,'numbertitle','off','name','level = 1')
% %imagesc(rec1(:,:,k));colorbar;
% mesh(x,y,rec1(:,:,k)); colorbar; 
% 
% norm_2 = norm(F - rec1(:,:,k),2);
% norm_1 = norm(F - rec1(:,:,k),1);
% 
% fprintf('The 2-norm difference is %0.8f.\n',norm_2)
% fprintf('The 1-norm difference is %0.8f.\n',norm_1)


























