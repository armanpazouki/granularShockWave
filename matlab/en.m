function [Etot, T, V] = en(vel)


V = zeros(501,1);
T = zeros(501,1);
Etot = zeros(501,1);

g = 9.81;

formatSpec = '%f,%f,%f,';
sizeA = [3 1];

for ii=1:501
    fileName = sprintf('./ex_vel_%03d/Energy/energy%04d.dat',vel,ii);

    fileID = fopen(fileName,'r');

    Kin = 0;
    Pot = 0;
    for jj = 1:101*100
       
        A = fscanf(fileID,formatSpec,sizeA);
        size(A);
        
        h = A(1,1);
        v = A(2,1);

        Pot = Pot + g*h;
        Kin = Kin + 0.7 * v^2;
    end
    
    V(ii,1) = Pot;
    T(ii,1) = Kin;
    Etot(ii,1) = Pot + Kin;
    
    fclose(fileID);
end

B = [Etot,V, T];
disp(B)

    


end