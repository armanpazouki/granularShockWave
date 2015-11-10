function zaleznosc(scal005_x, scal050_x, scal100_x, scal200_x,...
    scal005_z, scal050_z, scal100_z, scal200_z)

SavePlots = 1;
%switch 
% ii = 1: particle #50
% ii = 2: particle #60
% ii = 3: particle #70
% ii = 4: particle #80
% ii = 5: particle #90
ii = 1;
par = [50, 60, 70, 80, 90];
vel = [5, 50, 100, 200];

for ii=1:5

t_x = [scal005_x(ii,2), scal050_x(ii,2),...
    scal100_x(ii,2),scal200_x(ii,2)];

fig1 = figure(ii);
fig1.Position =[50 50 1280+50 720+50];

subplot(2,1,1)

plot(vel, t_x);
xlabel('Exciting velocity [m/s]');
ylabel('Time instant [s]')
plotTitle = sprintf('Time instant when the front of the wave reaches...\n... the particle #%d',par(ii)) ;
        
title(plotTitle) ;

%====================================

%switch 
% jj = 1: surface #50 
% jj = 2: surface #10
% jj = 3: surface #20
% jj = 4: surface #30
% jj = 5: surface #40
% jj = 6: surface #1 <- it is the same as par #1
jj = ii;
lay=[50, 10, 20, 30, 40, 1 ];

t_z = [scal005_z(ii,2), scal050_z(ii,2),...
    scal100_z(ii,2),scal200_z(ii,2)];

subplot(2,1,2)
plot(vel, t_z);
xlabel('Exciting velocity [m/s]');
ylabel('Time instant [s]')
plotTitle = sprintf('...the layer #%d',lay(ii));
title(plotTitle);

if (SavePlots == 1)
    figureName = sprintf('./time_vs_vel/time_vs_vel_part%04d',par(ii));%overlap
    saveas(fig1, figureName , 'png');
    saveas(fig1, figureName , 'fig');
end

end
end