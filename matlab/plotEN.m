function plotEN(Etot, vel)

    SavePlots = 1;

    fig1 = figure(1);
    fig1.Position =[50 50 1280+50 720+50];
    
    time = 1e-13*1:500;
    
    plot(time',Etot);
    xlabel('time [s]');
    ylabel('Total energy [J/kg]')
    plotTitle = sprintf('Total energy in the system. Exc. vel.: %d m/s',vel);
    title(plotTitle);
    
    if (SavePlots == 1)
        figureName = sprintf('./energy/energy_vel%03d',vel);%overlap
        saveas(fig1, figureName , 'png');
        saveas(fig1, figureName , 'fig');
    end


end