function F=RigidForce(X,Y)
global K;
F = K*(Y-X);
