function tor=torque(Z,ZCM,F)
global Nb

R = Z-repmat(ZCM,Nb,1);
tor = sum(R(:,2).* F(:,1) - R(:,1).* F(:,2));