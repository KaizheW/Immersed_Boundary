function F=ForceFilament(X,Z)
global Nb kp km ds K Kt;

T = K*(sqrt((X(kp,1)-X(km,1)).^2+(X(kp,2)-X(km,2)).^2)/(2*ds)-1);
T(1) = K*(sqrt((X(2,1)-X(1,1)).^2+(X(2,2)-X(1,2)).^2)/(ds)-1);
T(Nb) = K*(sqrt((X(Nb,1)-X(Nb-1,1)).^2+(X(Nb,2)-X(Nb-1,2)).^2)/(ds)-1);
tau = (X(kp,:)-X(km,:))./sqrt((X(kp,1)-X(km,1)).^2+(X(kp,2)-X(km,2)).^2);
tau(1,:) = (X(2,:)-X(1,:))./sqrt((X(2,1)-X(1,1)).^2+(X(2,2)-X(1,2)).^2);
tau(Nb,:) = (X(Nb,:)-X(Nb-1,:))./sqrt((X(Nb,1)-X(Nb-1,1)).^2+(X(Nb,2)-X(Nb-1,2)).^2);
Ttau = [T T].*tau;
F = (Ttau(kp,:)-Ttau(km,:))/(2*ds);
F(Nb,:) = (Ttau(Nb,:)-Ttau(Nb-1,:))/(ds);
F(1,:) = Kt*(Z-X(1,:));
