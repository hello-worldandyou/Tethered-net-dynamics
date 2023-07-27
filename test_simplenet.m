%缺点：时间过长就乱了、再长些成为一条线    运动时网格边长超过原边长很多       
%感觉步长不能取太大，因此不适合大型绳网 ，当然跟初速度大小也有关系

% Tension_calculation.m     计算绳网节点的受力和质量
% Tension_in_node.m         计算牵引绳对节点的拉力
% Tension_in_thread.m       计算绳网横线和竖线段的拉力

clear;
clc;
close all;
% kn1 kn2 是计算牵引绳的受力系数， k1 k2是计算绳网线段的受力系数 
global kn1 kn2 k1 k2 N dL0 mi mn

L = 5;                  % 绳网总长度
dL0 = 1;                % 绳网单元长度
Lct0 = 1.414*dL0;       % 牵引绳长度
N = L/dL0+1;            % 节点数  6×6
delta_t = 0.0001;       % 时间步长影响很大，步长太大，容易产生太大拉力？0.0002就网子就发散了，0.0001可以
tf = 1;                 % 40s       
Ns = round(tf/delta_t+1);      % 总步长
xnode = cell(Ns,1);     % cell数组一般被叫做元胞数组，它的每个单元可以储存不同的数据类型，可以是数值，字符或矩阵或元胞数组等，类似于学过的c语言里的结构体
xnode4 = cell(Ns,1);

rho_thread = 1440;                         % kg/m3
Enet = 7e10;                               % 杨氏模量  70GPA     
epsi= 0.106;                               % 阻尼比：damping ration of the net's thread
r_tether = 0.0001;                         % 绳子半径
r_ct = 0.002;                              % 牵引绳半径
mi = pi*r_tether^2*dL0*rho_thread;         % 绳段质量  
k1 = Enet*pi*r_tether^2/dL0;               % 弹力系数：k=N/m  Enet=N/m^2   EA/L0
k2 = 2*epsi*sqrt(mi*k1);                   % 阻尼系数：cij=2*epsi(zeta)*sqrt(m*k)

mn = pi*r_ct^2*Lct0*rho_thread;            % 牵引绳段的质量
kn1 = Enet*pi*r_ct^2/Lct0;
kn2 = 2*epsi*sqrt(mn*kn1);

m=0.5;                                      % 4个牵引航天器的质量
m_sate = m + mn/2;                          % update：+牵引绳半段的质量

T = zeros(3*N,N);
n0 = 0;                                     % 轨道角速度 sqrt(398600/42164^3)
Ao = [zeros(3,3),eye(3)
    3*n0^2,0,0,0,2*n0,0
    0,0,0,-2*n0,0,0
    0,0,-n0^2,0,0,0];                    % cw方程 
Bo = [zeros(3,3);eye(3)];                % 绳子节点方程  ，节点有3种质量
Bo1 = [zeros(3,3);eye(3)];              
dx = @(t,xi,u)(Ao*xi+Bo*u);              % u  代表 加速度
dx1 = @(t,xi,u)(Ao*xi+Bo1*u);



%% initial position of net node
xnode{1} = zeros(6*N,N);           % 通过二维复合矩阵表示三维矩阵挺好
xL0 = [10;0;L];                    % yz平面的第二象限 顶点   确定绳网左上角坐标初始位置，a点
x0 = cell(N,N);

for i = 1:N
    for j = 1:N
        xnode{1}((i-1)*6+1:(i-1)*6+3,j) = [xL0(1);xL0(2)+0.9*dL0*(j-1);xL0(3)-0.9*dL0*(i-1)];  %0.9表示未充分展开    % 初始位置确定，速度都是0(在轨道坐标系下也合理) 在一个平面，x轴坐标相同
%       xsave{i,j}(1:3,1) = [xL0(1);xL0(2)+dL0*(j-1);xL0(3)-dL0*(i-1)];
    end
end

% 拖曳航天器初始位置速度
xnode4{1}=zeros(6,4);

v=0;  
% 拖曳航天器 初始值   初始速度 1m/s
xnode4{1}(:,1)=[xL0(1)+Lct0;xL0(2)+dL0*(1-1);xL0(3)-dL0*(1-1); 1;-v;v];
xnode4{1}(:,2)=[xL0(1)+Lct0;xL0(2)+dL0*(N-1);xL0(3)-dL0*(1-1); 1;v;v];
xnode4{1}(:,3)=[xL0(1)+Lct0;xL0(2)+dL0*(N-1);xL0(3)-dL0*(N-1); 1;v;-v];
xnode4{1}(:,4)=[xL0(1)+Lct0;xL0(2)+dL0*(1-1);xL0(3)-dL0*(N-1); 1;-v;-v];
%% 积分
for k = 1:Ns-1
    
    T4 = Tension_in_node(xnode{k},xnode4{k});  %4个顶点与拖曳航天器之间的力
    [T1,T2] = Tension_in_thread(xnode{k});

    for i = 1:N
        for j = 1:N
            [Tij,mij] = Tension_calculation(i,j,T1,T2,T4);
            T(3*(i-1)+1:3*i,j) = Tij;
            m_net(i,j) = mij;
        end
    end
    
    for i = 1:N
        for j = 1:N
            xi = xnode{k}((i-1)*6+1:i*6,j);
            u = T(3*(i-1)+1:3*i,j)/m_net(i,j);                  %u 加速度啦
            [t,y] = ode45(@(t,xi)dx(t,xi,u),[(k-1)*delta_t,k*delta_t],xi);
            xi = y(end,:)';
            xnode{k+1}((i-1)*6+1:i*6,j) = xi;
%             xsave{i,j}(:,k+1) = xi;
        end
    end
   
 %% 4顶点  
 
 
 %xnode42(k,:)=xnode4(:,2)';    %拖曳航天器1时间历程图  可以放到下面循环里{}，把其他航天器也保存下来
 
 for i=1:4
        xi4 = xnode4{k}(:,i);
        u = -T4(:,i)/m_sate;
        [t,y] = ode45(@(t,xi)dx1(t,xi,u),[(k-1)*delta_t,k*delta_t],xi4);
        xi4 = y(end,:)';
        xnode4{k+1}(:,i) = xi4;       
 end 
        
    k
    
end
% 绘图
figure('color',[1 1 1]);        %图形背景改为白色；
for k=1:Ns
    cla;
    
        xx = zeros(N,N);
        yy = zeros(N,N);
        zz = zeros(N,N);
for i = 1:N
    for j = 1:N
        xx(i,j) = xnode{k}((i-1)*6+1,j);
        yy(i,j) = xnode{k}((i-1)*6+2,j);
        zz(i,j) = xnode{k}((i-1)*6+3,j);
    end
end

for i = 1:N
        plot3(xx(i,:),yy(i,:),zz(i,:),'b.',xx(i,:),yy(i,:),zz(i,:),'b');hold on;
        hold on;  
end

for j = 1:N
        plot3(xx(:,j),yy(:,j),zz(:,j),'b.',xx(:,j),yy(:,j),zz(:,j),'b');
        hold on;
end
        xnet=zeros(3,4);
        xnet(:,1)=xnode{k}(1:3,1);
        xnet(:,2)=xnode{k}(1:3,N);
        xnet(:,3)=xnode{k}((N-1)*6+1:(N-1)*6+3,N);
        xnet(:,4)=xnode{k}((N-1)*6+1:(N-1)*6+3,1);

for i=1:4
        plot3(xnode4{k}(1,i),xnode4{k}(2,i),xnode4{k}(3,i),'ro','MarkerFaceColor','r');
        hold on;
        x1=[xnet(1,i),xnode4{k}(1,i)];
        y1=[xnet(2,i),xnode4{k}(2,i)];
        z1=[xnet(3,i),xnode4{k}(3,i)];
        plot3(x1,y1,z1,'g')

end
        grid on;
        axis equal;
        axis([5 25 -5 10 -5 10])
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if k==1
         imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,'DelayTime',1e-8);
    else
         imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',1e-8);
    end

end

