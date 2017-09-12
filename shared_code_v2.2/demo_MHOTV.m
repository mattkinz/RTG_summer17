% demo to reproduce results from the article:
% T. Sanders, R.B. Platte,
% "Multiscale Higher Order TV Operators and Their Relationship to 
% Daubechies Wavelets"
%
% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 12/1/2016



clear
S_rates = .1:.1:0.7;  % sampling rates
mu = 200;   % regularization parameter
noise = false; % noise, true or false
d = 1024;   % dimension of the reconstruction
Ntrials = 20;   % Number of trials for each sampling rate
Nrates = numel(S_rates);    % number of sampling rates
cnt = 0;  % counter
tol = 1e-2;  % error tolerance to indicate exact recovery
njumps = 5; % number of jumps in the signal
porder = 2; % degree of the piecewise polynomial signal

%%
% store the recovery results in these matrices
Q = zeros(Ntrials,Nrates,3,'logical');
Qw = zeros(Ntrials,Nrates,3,'logical');
Qr = zeros(Ntrials,Nrates,3,'logical');

Q2 = zeros(Ntrials,Nrates,3);
Qw2 = zeros(Ntrials,Nrates,3);
Qr2 = zeros(Ntrials,Nrates,3);

% determine one random example to store and plot
[~,iii] = min(abs(S_rates-.2));
jjj = randi(Ntrials);
srec = zeros(d,4,3);  % store the results from the one example in this matrix

% loop for the simulations
for S_rate = S_rates
    cnt = cnt + 1;
for i = 1:Ntrials
    
    % generate the signal
    x = zeros(d,1);
    v = sort(randi(d,njumps,1));
    for p = 1:numel(v)-1
        pp = randi(2,1);
        c = randn(1)*3;
        b = randn(1)*3;
        a = randn(1)*3;
        if porder<1
            b = 0;
            a = 0;
        elseif porder <2
            a = 0;
        end
        alpha = randn(1)*3;
        beta = randn(1)*3;
        t = linspace(alpha,beta,v(p+1) - v(p) + 1);
        x(v(p):v(p+1)) = a*t.^2 + b*t + c;
    end
    nrmx = norm(x);
   

        A = sprand(round(d*S_rate),d,1/10);% random sampling matrix
        b = A*x;  % data vector
        if noise
           b = b + (n_levels(cnt2)*d/1024)*randn(d*S_rate,1); 
        end

        % options for l1 optimization algorithm
        clear pat;clear opt;
        pat.nonneg = false;
        pat.levels = 3;
        pat.inner_iter = 20;
        pat.outer_iter = 30;
        pat.tol_inn = 1e-5;
        pat.tol_out = 1e-5;
        pat.disp = false;
        pat.isreal = true;
        pat.data_mlp = true;
        
        % compute reconstructions for each order and each method (HOTV,
        % MHOTV and Daub wavelets)
        for k = 1:3     % loop over the order
            
            % HOTV reconstruction
            pat.order = k;
            pat.mu = mu;
            pat.levels = 1;
            rec1 = HOTV3D(A,b,[d,1,1],pat);
            if norm(rec1-x)/nrmx < tol
                Qr(i,cnt,k) = true;
            end
            Qr2(i,cnt,k) = norm(rec1-x)/nrmx;
            
            
            % MHOTV reconstruction, set pat.levels to 3
            pat.levels = 3;
            pat.mu = mu;
            rec2 = HOTV3D(A,b,[d,1,1],pat);
            if norm(rec2-x)/nrmx < tol
                Q(i,cnt,k) = true;
            end
            Q2(i,cnt,k) = norm(rec2-x)/nrmx;

            % Daubechies wavelets reconstruction
            wname = ['db',num2str(k)];
            opt = pat;
            [opt.D,opt.Dt] = my_wav_trans_1D(wname,pat.levels);
            
            rec3 = l1optimo(A,b,[d,1,1],opt);
            
            if norm(rec3-x)/nrmx < tol
                Qw(i,cnt,k) = true;
            end
            Qw2(i,cnt,k) = norm(rec3-x)/nrmx;
            
            
            % save one random example for demonstration
            if iii == cnt && jjj == i
                srec(:,1,k) = rec1;
                srec(:,2,k) = rec2;
                srec(:,3,k) = rec3;
                srec(:,4,k) = x;
            end
        end
      

        
end
end


%%
figure(3);hold off;
tnames = {'HOTV','MHOTV','db wavelets'};
for k = 1:3
    figure(3);subplot(2,2,k);hold off;
    plot(srec(:,4,1),'-','linewidth',1.5);hold on;
    plot(srec(:,k,1),'--','linewidth',1.5);
    plot(srec(:,k,3),'-.','linewidth',1.5);
    legend({'true','order 1 rec','order 3 rec'},'Location','best');
    title(tnames{k});
    xlabel('grid point');
    ylabel('function value');
    set(gca,'fontweight','bold','fontsize',12);
    axis([1 1024,min(srec(:,4,1)),max(srec(:,4,1))])
    hold off;
end

%%
% plot the results
for k = 1:3
    figure(1);subplot(3,1,k);hold off;
    plot(S_rates,sum(Qr(:,:,k),1)/Ntrials,'-.or','linewidth',1.5);hold on;
    plot(S_rates,sum(Q(:,:,k),1)/Ntrials,'--xb','linewidth',1.5);
    plot(S_rates,sum(Qw(:,:,k),1)/Ntrials,'-sk','linewidth',1.5);
    if k~=1
        legend({'HOTV','MHOTV','db'},'location','southeast');
    else
        legend({'TV','MHOTV','Haar'},'location','southeast');
    end
    xlabel('sampling rate');
    ylabel('probability of success');
    title({['Order ',num2str(k),' reconstruction models'];...
        ['Order ',num2str(porder),' piecewise polynomials']})
    grid on
    axis([min(S_rates),max(S_rates),0,1]);
    set(gca,'fontsize',10);
    %print(['constant-order',num2str(k)],'-dpng')
end


%%
for k = 1:3
    figure(2);subplot(3,1,k);hold off;
    plot(S_rates,sum(Qr2(:,:,k),1)/Ntrials,'-.or','linewidth',1.5);hold on;
    plot(S_rates,sum(Q2(:,:,k),1)/Ntrials,'--xb','linewidth',1.5);
    plot(S_rates,sum(Qw2(:,:,k),1)/Ntrials,'-sk','linewidth',1.5);
    if k~=1
        legend({'HOTV','MHOTV','db'},'location','northeast');
    else
        legend({'TV','MHOTV','Haar'},'location','northeast');
    end
    xlabel('sampling rate');
    ylabel('average error');
     title({['Order ',num2str(k),' reconstruction models'];...
         ['Order ',num2str(porder),' piecewise polynomials']})
    grid on
    axis([min(S_rates),max(S_rates),0,max(col(Qr2))]);
    set(gca,'fontsize',10);
    %print(['quadratic-mean-error-order',num2str(k)],'-dpng')
end
