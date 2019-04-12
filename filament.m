%% Immersed Boundary Method, 2D
% Flexible ribbon filament;

%% Initialize Parameter
clear
global Lx Ly Nx Ny Ks Kb Kt rho M mu g;
global h ipx ipy imx imy Nb ds kp km;
global a;
movie_or_not = 0; % whether export movie; 1->yes; 0->no.

% Global parameters
Lx = 1.0; % x size
Ly = 2.0; % y size
Nx = 64; % x mesh size
Ny = Nx/Lx*Ly; % y mesh size
Ks = 10.0; % stretch coefficient
Kb = 0.1; % bending rigidity
Kt = 500000; % target point pull back force constant
rho = 3*10^(-4); % fluid density
M = 4*10^(-4); % filament density
mu = 1*10^(-5); % fluid viscosity
g = 980; % gravity
tmax = 1; % time range
dt = 0.00001; % discretize time
clockmax = ceil(tmax/dt);

% Mesh
h = Lx/Nx; % grid size
ipx = [(2:Nx),1];
ipy = [(2:Ny),1];
imx = [Nx,(1:(Nx-1))];
imy = [Ny,(1:(Ny-1))];

% parameters specific for this code: filament.
L = 0.1*Ly; % length of the filament
Nb = ceil(L/(h/2)); 
ds = h/2;
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
ZX = Lx/2;
ZY = 3*Ly/4;
alpha = -pi/2; % initial tilted angle

% parameters specific for flow field.
u0 = 1.0; % initial uniform flow field velocity
dvorticity = 3; % delta vorticity, used to plot vorticity field
values= (-50*dvorticity):dvorticity:(50*dvorticity);

if movie_or_not == 1
    video = VideoWriter('targeted_filament.mp4','MPEG-4');
    video.FrameRate = 25;
    open(video);
end

%% Initialize Boundary and Flow Field
% generate a filament
X = zeros(Nb,2); % Boundary points
X(:,1) = ZX + ds*(0:(Nb-1))*cos(alpha);
X(:,2) = ZY + ds*(0:(Nb-1))*sin(alpha);
Z = X(1,:); % Fix the first point;
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
  u(1:2,:,1) = u0;
  u(1:2,:,2) = 0;
  XX=X+(dt/2)*interp(u,X);
  ff=spread_Filament(ForceFilament(XX,Z),XX);
  [u,uu]=fluid(u,ff);
  X=X+dt*interp(uu,XX);
  
  % Animation
  if mod(clock,1000)==1
      vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
      disp(max(max(vorticity))-min(min(vorticity)));
      if clock == 1
          dvorticity=(max(max(vorticity))-min(min(vorticity)));
          values= (-10000*dvorticity):20*dvorticity:(10000*dvorticity);
    %       valminmax=[min(values),max(values)];
      end
      contour(x,y,vorticity,values)
      colormap cool
      hold on
      plot(mod(XX(:,1),Lx),mod(XX(:,2),Lx),'ko')
      axis([0,Lx,0,Ly])
%       caxis(valminmax)
      axis equal
%       axis manual
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
