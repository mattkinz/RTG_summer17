function final_disp(out,opts)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

fprintf('Number of total iterations is %d. \n',out.total_iter);
fprintf('||Au-b||/||b||: %5.3f\n',out.final_error);
fprintf('Final ||w||_1: %5.3g\n',out.final_wl1);
fprintf('Final ||Du-w||^2: %5.3g\n',out.final_Du_w);
if opts.disp_fig
    if opts.disp
        figure(131);
        subplot(2,1,1);
        plot(out.lam1,'Linewidth',2);  %plot lam1, ||W||_1
        hold on;
        plot(out.lam3.*out.mu,'Linewidth',2);  %plot lam3, mu||Au -f||^2
        plot(abs(out.f),'Linewidth',2);   %plot f, the objective function
        plot(2:opts.inner_iter:max(size(out.f)),...
            out.f(2:opts.inner_iter:end),'kx','Linewidth',2);
        legend({'||W||_1','mu||Au - f||_2^2',...
            'obj function','mlp update'},...
            'Location','eastoutside');
        xlabel('iteration');
        grid on;
        set(gca,'fontweight','bold');
        hold off;

        subplot(2,1,2);
        plot(out.rel_error,'Linewidth',2);
        hold on
        plot(out.rel_chg_inn,'Linewidth',2);
        plot(out.rel_lam2,'Linewidth',2);
        plot(2:opts.inner_iter:max(size(out.f)),...
            0,'kx','Linewidth',2);
        legend({'Rel error','Rel chg','Rel lam2','mlp update'},'Location','eastoutside');
        xlabel('iteration');
        grid on;
        set(gca,'fontweight','bold');
        hold off

    end
end