% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 06/08/2016


% Running this script will reproduce some 1-D results related to the
% following publication:
% T. Sanders et. al., 
% "Recovering Fine Details from Under-Resolved Tomographic 
% Data using Higher Order Total Variation L1 regularization"



N_tests = 10;  % number of simulations
n=500;   %signal dimension
m = n*.25;  %sampling count


clear pat;

% set options for L1 optimization algorithm
pat.outer_iter =10;
pat.inner_iter = 15;
pat.mu = 100;
pat.beta = 32;
pat.disp = false;
pat.nonneg = true;
pat.max_c = false;
pat.data_mlp = false;


% length of the signal jumps
j_length = 10;

%build a basic nx1 signal
xp = zeros(n,1);
%xp(n/2+1:end) = 1;
xp(n/2-j_length/2-12:n/2 + j_length/2-12)=1;
xp(n/2-j_length/2+12:n/2 + j_length/2+12)=1;
%figure(1);hold off
%plot(xp(ind));hold on;
    

% initialize sampling matrix and data vectors
A = zeros(m,n,N_tests);
b = zeros(m,N_tests);
rec = zeros(n,3);


recs = cell(3,1);
% perform the reconstruction for each order of the finite difference
% operator
for k = 1:3
    pat.order = k;  % set the order




    %noise and sampling
    noise_mu = 0;   %mean noise level
    noise_sigma =.05;     %standard deviation


    recs{k} = zeros(n,N_tests);
    for ntests = 1:N_tests
        if k==1
            %Generate the random sampling matrix 
            B = rand(m,n);
            [~,s2] = max(B,[],2);
            %B = max(B-(1-sampling_perc/n),0);
            %A = zeros(m,n);
            for ii = 1:m
                %v = find(B(ii,:));
                v = s2(ii):min(s2(ii)+round(sqrt(n)),n);
                A(ii,v,ntests) = random('normal',1,.1,1,max(size(v)));
            end
            %generate sampling of x with A
            b(:,ntests) = A(:,:,ntests)*xp;
            eps = random('normal',noise_mu,noise_sigma,size(b,1),1);
            b(:,ntests) = b(:,ntests) + eps;
        end
        % run optimization 
        y = HOTV3D(A(:,:,ntests),b(:,ntests),[n,1,1],pat);


        recs{k}(:,ntests) = y;
    end
    rec(:,k) = mean(recs{k},2);
end



% plot the results
ind = n/2-50:n/2+50;
figure(1);hold off;
plot(ind,xp(ind));
hold on;

for k = 1:3
    figure(1);
    plot(ind,rec(ind,k));
end
figure(1);
axis([ind(1) ind(end) -.1 1.2]);
h =legend({'True signal','TV',...
    'HOTV Order 2','HOTV Order 3'},'Location','northeast');
xlabel('grid index');
ylabel('value');
title(['Mean reconstructed value from ',num2str(N_tests),' simulations'])
set(gca,'fontsize',12,'fontweight','bold','Position',[.12 .2 .76 1/2]);
set(h,'fontsize',10,'fontweight','normal');

% visualize one of the sampling matrices
figure(2);imagesc(A(:,:,end));
title({'Visual of typical sampling matrix';'used for these simulations'});