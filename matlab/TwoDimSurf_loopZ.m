function scal = TwoDimSurf_loopZ(vel)

clc;

GenePlots = 1;
SavePlots = 1;

scal = zeros(6,5);
aa_pomocnicze = [1,10,20,30,40,50];
aa = [50,959,1969,2979,3989,4999];


for ii = 1:6
    fileName = sprintf('./ex_vel_%03d/Particles/particle%04d.dat',vel,aa(ii));

    fileID = fopen(fileName,'r');
    formatSpec = '%f,%f,%f,%f,%f,%f,%f,';
    sizeA = [7 10000];

    A = fscanf(fileID,formatSpec,sizeA);
    size(A);
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
        fig1 = figure(aa_pomocnicze(ii));
        fig1.Position =[50 50 1280+50 720+50];

        subplot(2,1,1)
        [hAx] = plotyy(myTime,dispX,[myTime,myTime],[dispZ,sqrt(dispX.^2+dispZ.^2)]);
        legend('X', 'Z', 'Mag')
        xlabel('time [s]');
        ylabel(hAx(1),'dispX [\mum]') % left y-axis
        ylabel(hAx(2),'dispZ, dispMag [\mum]') % right y-axis
        plotTitle = sprintf('Interior. Exciting velocity: %d m/s.\n Displacement of the central particle from the layer #%d',vel, aa_pomocnicze(ii));
        title(plotTitle);
        
        axes(hAx(1))
        [timeX , maxDispX ] = ginput(1);
        axes(hAx(2))
        [timeZ , maxDispZ ] = ginput(1);

        subplot(2,1,2)
        [hAx] = plotyy(myTime,velX,[myTime,myTime],[velZ,sqrt(velX.^2+velZ.^2)]);
        legend('X', 'Z', 'Mag')
        xlabel('time [s]');
        ylabel(hAx(1),'velX [m/s]') % left y-axis
        ylabel(hAx(2),'velZ, velMag [m/s]') % right y-axis
        plotTitle = sprintf('Velocity of the central particle from the layer #%d',aa_pomocnicze(ii));
        title(plotTitle);

        if (SavePlots == 1)
            figureName = sprintf('./ex_vel_%03d/particleInterior%04d',vel,aa_pomocnicze(ii));%overlap
            saveas(fig1, figureName , 'png');
            saveas(fig1, figureName , 'fig');
        end
    end
    
    %[dispMaxX,dispIdxX]=max(abs(dispX));
    %[dispMaxZ,dispIdxZ]=max(abs(dispZ));
    
    %[maxVel,velIdx]=max(abs(velZ));
    
%     scal(ii,:) = [aa_pomocnicze(ii),dispIdx, myTime(dispIdx),maxDisp,...
%         velIdx, myTime(velIdx),maxVel];

   %  scal(ii,:) = [aa_pomocnicze(ii),dispIdxX, myTime(dispIdxX),dispMaxX,...
   %      dispIdxZ, myTime(dispIdxZ),dispMaxZ];
    
     scal(ii,:) = [aa_pomocnicze(ii), timeX, abs(maxDispX), timeZ, abs(maxDispZ)];
    
end

% ========================================================================
fig2 = figure(100);
fig2.Position =[50 50 1280+50 720+50];

partNum = scal(:,1);
dispTimeX = scal(:,2);
dispMaxX = scal(:,3);
dispTimeZ = scal(:,4);
dispMaxZ = scal(:,5);


% partNum = scal(:,1);
% idxDispX = scal(:,2);
% timeDispX = scal(:,3);
% maxDispX = scal(:,4);
% idxDispZ = scal(:,5);
% timeDispZ = scal(:,6);
% maxDispZ = scal(:,7);

subplot(2,1,1)
        [hAx] = plotyy(partNum, dispMaxX, partNum, dispMaxZ);
        legend('X', 'Z')
        xlabel('Layer''s #');
        ylabel(hAx(1),'dispX [\mum]') % left y-axis
        ylabel(hAx(2),'dispZ [\mum]') % right y-axis
        plotTitle = sprintf('Interior. Exciting velocity: %d m/s.\n Max. disp. of the central particles from the certain layers.',vel);
        title(plotTitle);

        subplot(2,1,2)
        plot(partNum, dispTimeX, partNum, dispTimeZ);
        xlabel('Layer''s #');
        ylabel('Time [s]'); % left y-axis
        plotTitle = sprintf('Time instant of the front of the wave reaching the layer.');
        title(plotTitle);

        if (SavePlots == 1)
            figureName = sprintf('./ex_vel_%03d/scalingInterior',vel);%overlap
            saveas(fig2, figureName , 'png');
            saveas(fig2, figureName , 'fig');
        end

end