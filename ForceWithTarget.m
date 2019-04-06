function F=ForceWithTarget(X,Z)
global kp km dtheta K Kt;

F=K*(X(kp,:)+X(km,:)-2*X)/(dtheta*dtheta);
F(1,:) = Kt*(Z-X(1,:));
