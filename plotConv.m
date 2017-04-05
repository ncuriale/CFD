function plotConv(normRes, bcFlag, iter)

if bcFlag == 1 % Use the following plot format if running subsonic nozzle case
   
    figure(3)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 iter]);
    box('on');
    hold('all');
    plot(log10(normRes),'k-','LineWidth',1.5);
    grid on
    xlabel('Iterations','interpreter','latex','FontSize',13);
    ylabel('log(Residual)','interpreter','latex','FontSize',13);
    legend({'Convergence'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Subsonic Nozzle','interpreter','latex','FontSize',12);

elseif bcFlag == 2 % Otherwise, use this format for the transonic nozzle case
    
    figure(3)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 iter]);
    box('on');
    hold('all');
    plot(log10(normRes),'k-','LineWidth',1.5);
    grid on
    xlabel('Iterations','interpreter','latex','FontSize',13);
    ylabel('log(Residual)','interpreter','latex','FontSize',13);
    legend({'Convergence'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Transonic Nozzle','interpreter','latex','FontSize',12);
    
end