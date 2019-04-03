%% Immersed Boundary Method, 2D
% This script is the main program.

%% Initialize Parameter
clear
global L N K rho mu tmax dt XC YC RC;
global a;
global h ip im Nb dtheta kp km;

L = 1.0; % size
N = 256; % mesh size
K = 10000000.0; % force constant
rho = 1.0; % fluid density
mu = 0.01; % fluid viscosity
tmax = 10; 
dt = 0.00001;
XC = L/4;
YC = L/2;
RC = L/32;
u0 = 5.0;
Amp = L/32;
T = 1;

h = L/N; % grid size
ip = [(2:N),1];
im = [N,(1:(N-1))];
Nb = ceil(pi*N); % generate a circle.
dtheta = 2*pi/Nb;
kp = [(2:Nb),1];
km = [Nb,(1:(Nb-1))];
clockmax = ceil(tmax/dt);

video = VideoWriter('cylinder.avi');
video.FrameRate = 25;
open(video);

%% Initialize Boundary and Flow Field
% generate a circle of ribbon
X = zeros(Nb,2);
Y = zeros(Nb,2);
X(:,1) = XC + RC*cos((1:Nb)*dtheta);
X(:,2) = YC + RC*sin((1:Nb)*dtheta);
Y(:,1) = XC + RC*cos((1:Nb)*dtheta);
Y(:,2) = YC + RC*sin((1:Nb)*dtheta);
% Coordinates, [0 h 2h ... L-h]
% Matrix index, (1 2 ... N)

% velocity of fluid flow
u=zeros(N,N,2);
[y,x] = meshgrid(0:h:L-h,0:h:L-h);
% u(:,:,1) = sin(pi*y/L);
u(:,:,1) = u0;

% vorticity: v_x - u_y; contour plot vorticity.
vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
% dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
% values= (-10*dvorticity):dvorticity:(10*dvorticity);
% valminmax=[min(values),max(values)];

set(gcf,'double','on');
% contour(x,y,vorticity,values);
contour(x,y,vorticity)
hold on
plot(X(:,1),X(:,2),'ko');
axis([0,L,0,L]);
% caxis(valminmax);
axis equal
axis manual
drawnow;
hold off

%% 4D matrix, fluid solver
a = zeros(N,N,2,2); % fluid solver
a(:,:,1,1) = ones(N,N);
a(:,:,2,2) = ones(N,N);

for m1=0:(N-1)
  for m2=0:(N-1)
    if~(((m1==0)||(m1==N/2))&&((m2==0)||(m2==N/2)))
      t=(2*pi/N)*[m1;m2];
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

for m1=0:(N-1)
  for m2=0:(N-1)
    t=(pi/N)*[m1;m2];
    s=sin(t);
    a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)...
                    /(1+(dt/2)*(mu/rho)*(4/(h*h))*(s'*s));
  end
end


%% Calculation
for clock=1:clockmax
  u(1:2,:,1) = u0;
  u(1:2,:,2) = 0;
%   Y = Oscillate(clock, T, Amp);
  Y(:,1) = XC + RC*cos((1:Nb)*dtheta);
  Y(:,2) = YC + RC*sin((1:Nb)*dtheta) + Amp*sin(2*pi/T*clock*dt);
  XX=X+(dt/2)*interp(u,X);
  ff=spread(RigidForce(XX,Y),XX);
  [u,uu]=fluid(u,ff);
  X=X+dt*interp(uu,XX);
  
  % Animation
  if mod(clock,1000)==1
      vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
      if clock == 1
          dvorticity=(max(max(vorticity))-min(min(vorticity)))/100;
          values= (-100*dvorticity):dvorticity:(100*dvorticity);
    %       valminmax=[min(values),max(values)];
      end
      contour(x,y,vorticity,values)
      hold on
      plot(X(:,1),X(:,2),'ko')
      axis([0,L,0,L])
    %   caxis(valminmax)
      axis equal
    %   axis manual
      title(['time = ',num2str(clock*dt)])
      drawnow
      hold off
      writeVideo(video,getframe(gcf));
      disp(clock*dt);
  end
end
%%
close(video);
