%% Immersed Boundary Method, 2D
% Flexible circular ribbon;

%% Initialize Parameter
clear
global Lx Ly Nx Ny K Kt rho mu tmax dt;
global h ipx ipy imx imy Nb dtheta kp km;
global a;
movie_or_not = 1; % whether export movie; 1->yes; 0->no.

% Global parameters
Lx = 2.0; % x size
Ly = 1.0; % y size
Nx = 256; % x mesh size
Ny = Nx/Lx*Ly; % y mesh size
K = 1.0; % force constant of the ribbon
Kt = 500000; % target point pull back force constant
rho = 1.0; % fluid density
mu = 0.01; % fluid viscosity
tmax = 5; % time range
dt = 0.00001; % discretize time

h = Lx/Nx; % grid size
ipx = [(2:Nx),1];
ipy = [(2:Ny),1];
imx = [Nx,(1:(Nx-1))];
imy = [Ny,(1:(Ny-1))];
RC = Ly/8; % radius of the cylinder;
Nb = ceil(2*pi*RC/(h/2)); %!!!!!!!!!!
dtheta = 2*pi/Nb; 
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
clockmax = ceil(tmax/dt);

% parameters specific for this code: cylinder.
CX = Lx/4;
CY = Ly/2;
u0 = 1.0; % prescribed inlet velocity;

if movie_or_not == 1
    video = VideoWriter('targeted_circle.avi');
    video.FrameRate = 25;
    open(video);
end

%% Initialize Boundary and Flow Field
% generate a circle of ribbon
X = zeros(Nb,2); % Boundary points
X(:,1) = CX + RC*cos((1:Nb)*dtheta + pi);
X(:,2) = CY + RC*sin((1:Nb)*dtheta);
Z = X(1,:); % Fix the first point;
% Coordinates, [0 h 2h ... L-h]
% Matrix index, (1 2 ... N)

% velocity of fluid flow
u=zeros(Nx,Ny,2);
[y,x] = meshgrid(0:h:Ly-h,0:h:Lx-h);
% u(:,:,2) = sin(2*pi*x/(Ly));
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
plot(X(:,1),X(:,2),'ko');
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
a = zeros(Nx,Ny,2,2); % fluid solver
a(:,:,1,1) = ones(Nx,Ny);
a(:,:,2,2) = ones(Nx,Ny);

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
  ff=spread(ForceWithTarget(XX,Z),XX);
  [u,uu]=fluid(u,ff);
  X=X+dt*interp(uu,XX);
  
  % Animation
  if mod(clock,1000)==1
      vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
%       disp(max(max(vorticity))-min(min(vorticity)));
      if clock == 1
          dvorticity=(max(max(vorticity))-min(min(vorticity)));
          values= (-500*dvorticity):10*dvorticity:(500*dvorticity);
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
