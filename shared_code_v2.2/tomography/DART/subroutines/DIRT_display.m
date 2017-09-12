function [] = DIRT_display(U,opt,out,i,k)


    if i==1
        figure(1);
        subplot(2,2,2);imagesc(U,[opt.grays(1),opt.grays(end)]);title('Initial Segmentation');
        figure(1);
    else
        figure(1);
        subplot(2,2,3);imagesc(U,[opt.grays(1),opt.grays(end)]);
        figure(1);
        title(sprintf('DART iteration #%i',i));
    end
    figure(1);
    subplot(2,2,4);
    plot(1:opt.total_iter,...
        [out.iterative_error(1:i,k);zeros(opt.total_iter-i,1)],'-*');
    hold on
   % plot(1:opt.total_iter,...
    %    [out.set(2:i+1,k)/numel(U);zeros(opt.total_iter-i,1)],'g');
    %legend('Box','on',{'||Wx-b||/||b||','relative # of free pixels'});
    hold off
    
    figure(1);      
    if i==1
    title('Residual at each DART iteration');
    end
    
end