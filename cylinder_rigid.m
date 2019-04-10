%% Immersed Boundary Method, 2D
% Penalty Immersed Boundary Method, Rigid Body.

%% Initialize Parameter
clear
global Lx Ly Nx Ny K rho m mu g tmax dt;
global h ds ipx ipy imx imy Nb dtheta kp km;
global a;
movie_or_not = 0; % whether export movie; 1->yes; 0->no.

% Global parameters
Lx = 1.0; % x size
Ly = 2.0; % y size
Nx = 128; % x mesh size
Ny = Nx/Lx*Ly; % y mesh size
K = 500000.0; % force constant
rho = 1.0; % fluid density
m = 0.03; % Cylinder excess mass per unit length
mu = 0.01; % fluid viscosity
g = 980; % gravity
tmax = 1; % time range
dt = 0.00001; % discretize time
clockmax = ceil(tmax/dt);

% Mesh
h = Lx/Nx; % grid size
ds = h/2;
ipx = [(2:Nx),1];
ipy = [(2:Ny),1];
imx = [Nx,(1:(Nx-1))];
imy = [Ny,(1:(Ny-1))];

% parameters specific for this code: cylinder
RC = Lx/32; % radius of the cylinder;
Nb = ceil(2*pi*RC/(h/2)); % points at boundary
dtheta = 2*pi/Nb; % polar
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
ZCM = [Lx/2 7*Ly/8]; % Massive component initial mass center
VCM = [0 0]; % Massive component initial velocity
Omega = 0; % Massive component initial angular velocity

% if ellipse, more parameters: x^2/a^2 + y^2/b^2 = 1
Ratio = 2; % b/a
TiltAngle = pi/3; % initial tilted angle

% parameters specific for flow field.
u0 = 0.0; % initial uniform flow field velocity
dvorticity=5; % delta vorticity, used to plot vorticity field
values= (-30*dvorticity):dvorticity:(30*dvorticity);

% Movie export initial.
if movie_or_not == 1
    video = VideoWriter('cylinder_fall2.mp4','MPEG-4');
    video.FrameRate = 25;
    open(video);
end

%% Initialize Boundary and Flow Field
% generate a circle of ribbon
X = zeros(Nb,2); % Boundary points
% Z = zeros(Nb,2); % Target points
C = zeros(Nb,2); % Coordinates
X(:,1) = ZCM(1) + RC*cos((1:Nb)*dtheta+TiltAngle);
X(:,2) = ZCM(2) + Ratio*RC*sin((1:Nb)*dtheta);
Z = X;
C(:,1) = RC*cos((1:Nb)*dtheta+TiltAngle); % Z - ZCM
C(:,2) = Ratio*RC*sin((1:Nb)*dtheta);
I0 = m*(sum(sum(C.^2)))*ds;
M = Nb*ds*m;
L = I0*Omega;
E = eye(2);
% Coordinates, [0 h 2h ... L-h]
% Matrix index, (1 2 ... N)

% velocity of fluid flow
u=zeros(Nx,Ny,2);
[y,x] = meshgrid(0:h:Ly-h,0:h:Lx-h);
% u(:,:,1) = sin(pi*y/L);
u(:,:,1) = u0;

% vorticity: v_x - u_y; contour plot vorticity.
vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
% dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
% values= (-10*dvorticity):dvorticity:(10*dvorticity);
% valminmax=[min(values),max(values)];

figure('Position', [1 1 round(1000*Lx) round(1000*Ly)])
% set(gcf,'double','on');
% contour(x,y,vorticity,values);
contour(x,y,vorticity)
hold on
plot(X(:,1),X(:,2),'k.');
axis([0,Lx,0,Ly]);
% caxis(valminmax);
axis equal
axis manual
drawnow;
if movie_or_not == 1
    writeVideo(video,getframe(gcf));
end
hold off

%% 4D matrix, fluid solver
a = zeros(Nx,Ny,2,2); a(:,:,1,1) = ones(Nx,Ny); a(:,:,2,2) = ones(Nx,Ny);
for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    if~(((m1==0)||(m1==Nx/2))&&((m2==0)||(m2==Ny/2)))
      t=[2*pi/Nx;2*pi/Ny].*[m1;m2];
      s=sin(t);
      ss=(s*s')/(s'*s);
      a(m1+1,m2+1,1,1)=a(m1+1,m2+1,1,1)-ss(1,1);
      a(m1+1,m2+1,1,2)=a(m1+1,m2+1,1,2)-ss(1,2);
      a(m1+1,m2+1,2,1)=a(m1+1,m2+1,2,1)-ss(2,1);
      a(m1+1,m2+1,2,2)=a(m1+1,m2+1,2,2)-ss(2,2);
    end
  end
end
for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    t=[pi/Nx;pi/Ny].*[m1;m2];
    s=sin(t);
    a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)...
                    /(1+(dt/2)*(mu/rho)*(4/(h*h))*(s'*s));
  end
end

%% Calculation
for clock=1:clockmax
  XX=X+(dt/2)*interp(u,X);
  ZZCM = ZCM + (dt/2)*VCM;
  Omega = L/I0;
  theta = Omega*dt/2;
  RM = [cos(theta) -sin(theta);sin(theta) cos(theta)];
  EE = RM * E;
  ZZ = ZZCM + (EE*C')';
  FF = RigidForce(XX,ZZ);
  ff = spread(FF,XX);
  TT = torque(ZZ,ZZCM,FF);
  [u,uu]=fluid(u,ff);
  VVCM = VCM + dt/(2*M)*sum(-FF) - [0 dt*g/2];
  LL = L + dt*TT/2;
  X = X + dt*interp(uu,XX);
  ZCM = ZCM + dt*VVCM;
  OOmega = LL/I0;
  theta = OOmega*dt;
  RM = [cos(theta) -sin(theta);sin(theta) cos(theta)];
  E = RM * E;
  Z = ZCM + (E*C')';
  VCM = VCM + dt*sum(-FF)/M - [0 dt*g];
  L = L + dt*TT;
  
  % Animation
  if mod(clock,1000)==1
      vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
%       disp(Omega);
%       disp(max(max(vorticity))-min(min(vorticity)));
      contour(x,y,vorticity,values)
%       contour(x,y,vorticity)
      colormap cool
      hold on
      plot(mod(X(:,1),Lx),mod(X(:,2),Ly),'k.')
      hold on
      plot(mod(X(1,1),Lx),mod(X(1,2),Ly),'ro')
      axis([0,Lx,0,Ly])
    %   caxis(valminmax)
      axis equal
    %   axis manual
      title(['time = ',num2str(clock*dt)])
      drawnow
      hold off
      if movie_or_not == 1
        writeVideo(video,getframe(gcf));
      end
      disp(clock*dt);
  end
end
%%
if movie_or_not == 1
    close(video);
end
