function plotSol( P_d, M_d, rho_d,P_exact, M_exact, RHO_exact, bcFlag, nodes)

if (bcFlag == 1)
  
    %Mach
    figure(1)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(M_exact(2:nodes),'r-','LineWidth',1.5);
    plot(M_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Mach Number (-)','interpreter','latex','FontSize',13);
    legend({'Exact','Calculated'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Subsonic Nozzle','interpreter','latex','FontSize',12);
       
    %Pressure
    figure(2)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(P_exact(2:nodes),'r-','LineWidth',1.5);
    plot(P_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Pressure (Pa)','interpreter','latex','FontSize',13);
    legend({'Exact','Calculated'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Subsonic Nozzle','interpreter','latex','FontSize',12);

elseif (bcFlag == 2) 

    %Mach
    figure(1)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(M_exact(2:nodes),'r-','LineWidth',1.5);
    plot(M_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Mach Number (-)','interpreter','latex','FontSize',13);
    legend({'Exact','Calculated'},'interpreter','latex','FontSize',11,'Location','NorthWest')
    title('Transonic Nozzle','interpreter','latex','FontSize',12);
        
    %Pressure
    figure(2)
    axes('YMinorTick','on','YGrid','on','XMinorTick','on','XGrid','on');
    xlim([0 nodes]);
    box('on');
    hold('all');
    plot(P_exact(2:nodes),'r-','LineWidth',1.5);
    plot(P_d(1:nodes-1),'kx');
    grid on
    xlabel('Nodes','interpreter','latex','FontSize',13);
    ylabel('Pressure (Pa)','interpreter','latex','FontSize',13);
    legend({'Exact','Calculated'},'interpreter','latex','FontSize',11,'Location','NorthEast')
    title('Transonic Nozzle','interpreter','latex','FontSize',12);

end
    

