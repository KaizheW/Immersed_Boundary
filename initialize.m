global L N tmax dt;
global h ip im Nb dtheta;

% generate a circle of ribbon
clockmax = ceil(tmax/dt);
X = zeros(Nb,2);
X(:,1) = L/2 + L/4*cos((1:Nb)*dtheta);
X(:,2) = L/2 + L/4*sin((1:Nb)*dtheta);

% velocity of fluid flow
u=zeros(N,N,2);
[y,x] = meshgrid(0:h:L-h,0:h:L-h);
u(:,:,2) = sin(2*pi*x/L);

% vorticity: v_x - u_y;
vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
values= (-10*dvorticity):dvorticity:(10*dvorticity);
valminmax=[min(values),max(values)];

set(gcf,'double','on')
contour(x,y,vorticity,values)
hold on
plot(X(:,1),X(:,2),'ko')
axis([0,L,0,L])
caxis(valminmax)
axis equal
axis manual
drawnow
hold off

