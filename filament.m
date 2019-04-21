%% Immersed Boundary Method, 2D
% Flexible ribbon filament;

% Some ideas: Ktt should be greater than Kt, for at least an order;
% Bending regidity is important. If zero, huge vibration. Best value not
% find yet.

%% Initialize Parameter
clear
global Lx Ly Nx Ny Ks Kb Kt rho M mu g dt;
global h ipx ipy imx imy Nb ds kp km;
global a;
movie_or_not = 1; % whether export movie; 1->yes; 0->no.

% Global parameters
Lx = 10; % x size
Ly = 20; % y size
Nx = 128; % x mesh size
Ny = Nx/Lx*Ly; % y mesh size
Ks = 5e6; % stretch coefficient
Kb = 1e2; % bending rigidity
Kt = 5e6; % between massless and massive filament
% Ktt = 1e8; % target point
rho = 1; % fluid density
M = 1; % filament density
mu = 1e-2; % fluid viscosity
g = 0; % gravity
tmax = 1; % time range
dt = 1e-6; % discretize time
clockmax = ceil(tmax/dt);
alpha0 = 1000;

% Mesh
h = Lx/Nx; % grid size
ipx = [(2:Nx),1];
ipy = [(2:Ny),1];
imx = [Nx,(1:(Nx-1))];
imy = [Ny,(1:(Ny-1))];

% parameters specific for this code: filament.
L = 2; % length of the filament
Nb = ceil(L/(h/2))+1; 
ds = h/2;
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
ZX = Lx/2; % fixed point; first point of the filament
ZY = 13*Ly/16; % Y
alpha = -pi/2+0.1; % initial tilted angle; -pi/2 -> vertical down

% parameters specific for flow field.
u0 = -10.0; % initial uniform flow field velocity
dvorticity = 20; % delta vorticity, used to plot vorticity field
values= [(-100*dvorticity):dvorticity:(-1*dvorticity), ...
    (0*dvorticity):dvorticity:(100*dvorticity)];

if movie_or_not == 1
    video = VideoWriter('targeted_filament4.mp4','MPEG-4');
    video.FrameRate = 30;
    open(video);
end

%% Initialize Boundary and Flow Field
% generate a filament
X = zeros(Nb,2); % Boundary points
X(:,1) = ZX + ds*(0:(Nb-1))*cos(alpha);
X(:,2) = ZY + ds*(0:(Nb-1))*sin(alpha);
Y = X; % Massive boundary
Z = [ZX ZY]; % Fix the first point;
V = zeros(Nb,2);
% Coordinates, [0 h 2h ... L-h]
% Matrix index, (1 2 ... N)

% velocity of fluid flow
u=zeros(Nx,Ny,2);
[y,x] = meshgrid(0:h:Ly-h,0:h:Lx-h);
% u(:,:,2) = sin(2*pi*x/(Ly));
u(:,:,2) = u0;

% vorticity: v_x - u_y; contour plot vorticity.
vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
% dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
% values= (-10*dvorticity):dvorticity:(10*dvorticity);
% valminmax=[min(values),max(values)];

figure('Position', [1 1 round(100*Lx) round(100*Ly)])
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
%     a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)-(s*s')/(s'*s);
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
  u(:,end-1:end,1) = 0;
  u(:,end-1:end,2) = u0;
  XX = X + (dt/2)*interp(u,X);
  YY = Y + (dt/2)*V;
  FF1 = ForceFilament(XX,YY,Z);
  ff = spread_Filament(FF1,XX);
%   ff(:,end-1:end,1) = ff(:,end-1:end,1) + alpha0*(-u(:,end-1:end,1));
%   ff(:,end-1:end,2) = ff(:,end-1:end,2) + alpha0*(u0-u(:,end-1:end,1));
  [u,uu] = fluid(u,ff);
  FF2 = Kt*(YY-XX);
  VV = V + (-FF2)*(dt/2)/M;
  X = X + dt*interp(uu,XX);
  Y = Y + dt*VV;
  V = V + (-FF2)*dt/M;
  
  % Animation
  if mod(clock,1000)==0
      vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
      disp(max(max(vorticity))-min(min(vorticity)));
      contour(x,y,vorticity,values)
      colormap cool
      hold on
      plot(mod(X(:,1),Lx),mod(X(:,2),Ly),'k.')
      hold on
      plot(ZX, ZY,'ro')
      axis([0,Lx,0,Ly])
%       caxis(valminmax)
      axis equal
      axis manual
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
