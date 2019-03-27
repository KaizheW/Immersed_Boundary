%% Immersed Boundary Method, 2D
% This script is the main program.

%%
clear
global L N K rho mu tmax dt;
global a;
global h ip im Nb dtheta kp km;

L = 1.0; % size
N = 64; % mesh size
K = 1.0; % force constant
rho = 1.0; % fluid density
mu = 0.01; % fluid viscosity
tmax = 1; 
dt = 0.01;

h = L/N; % grid size
ip = [(2:N),1];
im = [N,(1:(N-1))];
Nb = ceil(pi*N); % generate a circle.
dtheta = 2*pi/Nb;
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
initialize;

a = zeros(N,N,2,2); % fluid solver
init_a;

%%
for clock=1:clockmax
  XX=X+(dt/2)*interp(u,X);
  ff=spread(Force(XX),XX);
  [u,uu]=fluid(u,ff);
  X=X+dt*interp(uu,XX);
  
  % Animation
  vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
  contour(x,y,vorticity,values)
  hold on
  plot(X(:,1),X(:,2),'ko')
  axis([0,L,0,L])
  caxis(valminmax)
  axis equal
  axis manual
  drawnow
  hold off
end
%%
