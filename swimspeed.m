%% Calculate the swimming speed of the magnetic swimmer
% Compare one swimmer and double swimmer

%% Control panel
p=[ 1;    % May14t1 and May6t1, Af=0.05, T=0.05, L=0.2
    0;    % May14t2 and May6t2, Af=0.05, T=0.1,  L=0.2
    0;    % May14t3 and May6t3, Af=0.05, T=0.05, L=0.3
    0;    % May14t4 and May6t4, Af=0.05, T=0.1,  L=0.3
    0;    % May14t5 and May6t5, Af=0.05, T=0.05, L=0.4
    0];   % May14t6 and May6t6, Af=0.05, T=0.1,  L=0.4

%% May14t1 (one) and May6t1 (two), Af=0.05, T=0.05, L=0.2
if p(1) == 1
    date = 'May14';
    trial = 1;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 1;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end

%% May14t2 (one) and May6t2 (two), Af=0.05, T=0.1, L=0.2
if p(2) == 1
    date = 'May14';
    trial = 2;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 2;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end

%% May14t3 (one) and May6t3 (two), Af=0.05, T=0.05, L=0.3
if p(3) == 1
    date = 'May14';
    trial = 3;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 3;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end

%% May14t4 (one) and May6t4 (two), Af=0.05, T=0.1, L=0.3
if p(4) == 1
    date = 'May14';
    trial = 4;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 4;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end

%% May14t5 (one) and May6t5 (two), Af=0.05, T=0.05, L=0.4
if p(5) == 1
    date = 'May14';
    trial = 5;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 5;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end

%% May14t6 (one) and May6t6 (two), Af=0.05, T=0.1, L=0.4
if p(6) == 1
    date = 'May14';
    trial = 6;
    maxstep = 1000;
    front0 = zeros(maxstep,2);
%     rear0 = zeros(maxstep,2);
    speed0 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X','ZX','ZY');
        front0(i,:) = X(1,:);
        speed0(i,:) = (X(1,:) - [ZX ZY])/time(i);
    end
    plot(time, speed0(:,2),'k','DisplayName','Single');
    hold on

    date = 'May6';
    trial = 6;
    maxstep = 1000;
    front1 = zeros(maxstep,2);
    % rear1 = zeros(maxstep,2);
    speed1 = zeros(maxstep,2);
    front2 = zeros(maxstep,2);
    % rear2 = zeros(maxstep,2);
    speed2 = zeros(maxstep,2);
    load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_1.mat'],'Af','T','dt');
    dt = dt*1000;
    time = (dt:dt:maxstep*dt)';

    for i=1:maxstep
    %     disp(i)
        load(['/Users/kaizhe/Desktop/',date,'t',num2str(trial),'/',date,'t',...
            num2str(trial),'_',num2str(i),'.mat'],'X1','X2','ZX1','ZX2','ZY1','ZY2');
        front1(i,:) = X1(1,:);
        front2(i,:) = X2(1,:);
        speed1(i,:) = (X1(1,:)-[ZX1 ZY1])/time(i);
        speed2(i,:) = (X2(1,:)-[ZX2 ZY2])/time(i);
    end

    plot(time, speed1(:,2), 'r','DisplayName','Double 1');
    hold on
    plot(time, speed2(:,2), 'b','DisplayName','Double 2');
    legend('show')
end