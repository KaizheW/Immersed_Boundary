function F=ForceFilament(X,Y)
global Nb ds Ks Kb Kt;

% s = 1, 2, 3, ..., Nb-1, Nb;
% size of T, tau and Ttau: Nb-1; Size of Fs, Fb, and F: Nb;
T = Ks*( sqrt(sum((X(2:Nb,:) - X(1:Nb-1,:)).^2, 2))/ds - 1 );
tau = (X(2:Nb,:) - X(1:Nb-1,:))./sqrt(sum((X(2:Nb,:) - X(1:Nb-1,:)).^2, 2));
Ttau = [T T].*tau;
Fs = ([Ttau; 0, 0] - [0, 0; Ttau])/ds;

Fb = zeros(Nb, 2);
Fb(1,:) = -X(1,:) + 2*X(2,:) - X(3,:);
Fb(2,:) = 2*X(1,:) - 5*X(2,:) + 4*X(3,:) - X(4,:);
Fb(3:Nb-2,:) = -X(1:Nb-4,:) + 4*X(2:Nb-3,:) - 6*X(3:Nb-2,:) + 4*X(4:Nb-1,:) - X(5:Nb,:);
Fb(Nb-1,:) = -X(Nb-3,:) + 4*X(Nb-2,:) - 5*X(Nb-1,:) + 2*X(Nb,:);
Fb(Nb,:) = -X(Nb-2,:) + 2*X(Nb-1,:) - X(Nb,:);
Fb = Kb*Fb/(ds^4);

Ft = Kt*(Y-X);

F = Fs + Fb + Ft;
% F(1,:) = F(1,:) + Kt*(Z-X(1,:));
