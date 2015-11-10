function scal = TwoDimSurf_loopX(vel)

clc;

GenePlots = 1;
SavePlots = 1;

%read time
% myTimeFileName = sprintf('../time.dat');
% myTimeFileID = fopen(myTimeFileName,'r');
% myTimeFormatSpec = '%f,';
% myTimeSize = [1 Inf];
% myTime = fscanf(myTimeFileID,myTimeFormatSpec,myTimeSize);

% read particle 51

scal = zeros(5,5);

%aa = [50,60,70,80,90,100];% to 100th particle wave does not come
aa = [50,60,70,80,90];


for ii = 1:5
    fileName = sprintf('./ex_vel_%03d/Particles/particle%04d.dat',vel,aa(ii));

    fileID = fopen(fileName,'r');
    formatSpec = '%f,%f,%f,%f,%f,%f,%f,';
    sizeA = [7 10000];

    A = fscanf(fileID,formatSpec,sizeA);
    size(A) ;
    fclose(fileID);

    myTime = A(1,:)';
    posX = A(2,:)';
    posZ = A(3,:)';
    angY = A(4,:)';
    dispX = posX-posX(1,1);
    dispZ = posZ-posZ(1,1);
    dangY = angY-angY(1,1);
    velX = A(5,:)';
    velZ = A(6,:)';
    avelY = A(7,:)';
    
    if (GenePlots == 1)
        fig1 = figure(aa(ii));
        fig1.Position =[50 50 1280+50 720+50];

        subplot(2,1,1)
        [hAx] = plotyy(myTime,dispX,[myTime,myTime],[dispZ,sqrt(dispX.^2+dispZ.^2)]);
        if (ii>1)
            legend('X', 'Z', 'Mag','Location','NorthWest')
        else
            legend('X', 'Z', 'Mag','Location','NorthEast')
        end
        xlabel('time [s]');
        ylabel(hAx(1),'dispX [\mum]') % left y-axis
        ylabel(hAx(2),'dispZ, dispMag [\mum]') % right y-axis
        plotTitle = sprintf('Surface. Exciting velocity: %d m/s.\nDisplacement of the particle #%d',vel, aa(ii));
        title(plotTitle);
        
        axes(hAx(1))
        [timeX , maxDispX ] = ginput(1);
        axes(hAx(2))
        [timeZ , maxDispZ ] = ginput(1);

        subplot(2,1,2)
        [hAx] = plotyy(myTime,velX,[myTime,myTime],[velZ,sqrt(velX.^2+velZ.^2)]);
        if (ii>1)
            legend('X', 'Z', 'Mag','Location','NorthWest')
        else
            legend('X', 'Z', 'Mag','Location','NorthEast')
        end
        xlabel('time [s]');
        ylabel(hAx(1),'velX [m/s]') % left y-axis
        ylabel(hAx(2),'velZ, velMag [m/s]') % right y-axis
        plotTitle = sprintf('Velocity of the particle #%d',aa(ii));
        title(plotTitle);

        if (SavePlots == 1)
            figureName = sprintf('./ex_vel_%03d/particleSurface%04d',vel,aa(ii));%overlap
            saveas(fig1, figureName , 'png');
            saveas(fig1, figureName , 'fig');
        end
    end
    
    %[maxDispX,dispIdxX]=max(abs(dispX)); % getting that by clicking on the graph
    %[maxDispZ,dispIdxZ]=max(abs(dispZ)); % getting that by clicking on the graph
    %[maxVel,velIdx]=max(abs(velX)); % obsolete
    
%    scal(ii,:) = [aa(ii),dispIdx, myTime(dispIdx),maxDisp,...
%        velIdx, myTime(velIdx),maxVel];
 
%     scal(ii,:) = [aa(ii),dispIdxX, myTime(dispIdxX),maxDispX, ...
%         dispIdxZ, myTime(dispIdxZ),maxDispZ];

    scal(ii,:) = [aa(ii), timeX, abs(maxDispX), timeZ, abs(maxDispZ)];
    
end



% ========================================================================
fig2 = figure(1);
fig2.Position =[50 50 1280+50 720+50];

partNum = scal(:,1);
dispTimeX = scal(:,2);
dispMaxX = scal(:,3);
dispTimeZ = scal(:,4);
dispMaxZ = scal(:,5);



%dispTime - velTime % just to check - no sense you idiot!!!
% the wave front can be seen only at displacement plot!!!

subplot(2,1,1)
        [hAx] = plotyy(partNum, dispMaxX,partNum, dispMaxZ);
        legend('dispMaxX','dispMaxZ')
        xlabel('Particle''s #');
        ylabel(hAx(1),'dispX [\mum]') % left y-axis
        ylabel(hAx(2),'dispZ [\mum]') % left y-axis
        plotTitle = sprintf('Particles on the surface. Exciting velocity: %d m/s.\n MaxDispX, MaxDispZ.',vel);
        title(plotTitle);

        subplot(2,1,2)
        plot(partNum, dispTimeX, partNum, dispTimeZ);
        xlabel('Particle''s #');
        ylabel('Time [s]') % left y-axis
        legend('X','Z')
        plotTitle = sprintf('Time instant of dispMaxX, dispMaxZ');
        title(plotTitle);
        
        if (SavePlots == 1)
            figureName = sprintf('./ex_vel_%03d/scalingSurface',vel);%overlap
            saveas(fig2, figureName , 'png');
            saveas(fig2, figureName , 'fig');
        end

end