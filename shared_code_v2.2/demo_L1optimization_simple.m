% demo to reproduce results from the article:
% T. Sanders, R.B. Platte,
% "Multiscale Higher Order TV Operators and Their Relationship to 
% Daubechies Wavelets"
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 02/22/2017



d = 2048;  % signal dimension
S_rate = .5; % sampling rate
mu = 100;  % regularization parameter
noise = true;  % add noise to data?

%% build the simple piecewise polynomial signal
x = zeros(d,1);
x(round(d/6):round(d/6+d/4)) = ...
    2-linspace(-sqrt(2),sqrt(2),numel(round(d/6):round(d/6+d/4))).^2;
x(round(d/2):round(d/2+d/6)) = 1;
x(round(d/2+d/6):round(d/2+d/3)) = linspace(1,0,numel(round(d/2+d/6):round(d/2+d/3)));
%% generate random sampling matrix
A = sprand(d*S_rate,d,1/10);%randn(d*S_rate,d);
b = A*x;
if noise
   b = b + (0.5*d/1024)*randn(d*S_rate,1); 
end


%% set up l1 optimization parameters
clear pat;
pat.mu = mu;
pat.nonneg = false;
pat.levels = 1;
pat.inner_iter = 10;
pat.outer_iter = 50;
pat.tol_inn = 1e-7;
pat.tol_out = 1e-7;
pat.disp = false;
pat.isreal = true;
pat.data_mlp = false;



%% run the software to generate the results for multiple orders and multiple
% scales (MHOTV)
rec = zeros(d,3,3);
for levels = 1:3
    pat.levels = levels;
    for k = 1:3
        pat.order = k;
        [rec(:,k,levels),out] = HOTV3D(A,b,[d,1,1],pat);
    end
end


rcgls = run_cgs(A,b,1e-7,100);
%% display

% plot the regularized solutions
figure(2);hold off;
subplot(4,1,1);
plot(x,'-k','linewidth',2);
legend('true');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold');
%axis('off');

subplot(4,1,2);
plot(rec(:,1,1),'-k','linewidth',2);
legend('TV');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');

subplot(4,1,3);
plot(rec(:,2,1),'-k','linewidth',2);
legend('HOTV2');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');


subplot(4,1,4);
plot(rec(:,3,1),'-k','linewidth',2);
legend('HOTV3');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');



figure(3);hold off;
subplot(4,1,1);
plot(x,'-k','linewidth',2);
legend('true');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold');
%axis('off');

subplot(4,1,2);
plot(rec(:,1,3),'-k','linewidth',2);
legend('Haar');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');

subplot(4,1,3);
plot(rec(:,2,3),'-k','linewidth',2);
legend('MHOTV2');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');


subplot(4,1,4);
plot(rec(:,3,3),'-k','linewidth',2);
legend('MHOTV3');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);
%axis('off');


figure(4);hold off;
subplot(4,1,1);
plot(x,'-k','linewidth',2);
legend('true');
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold');
%axis('off');

subplot(4,1,2);
plot(rcgls,'-k','linewidth',2);
legend({'least squares'});
axis([0 d -0.5 2.5])
set(gca,'fontsize',12,'fontweight','bold','XTick',[],'YTick',[]);%,'Position',[.12 .2 .76 .5]);

