%ȱ�㣺ʱ����������ˡ��ٳ�Щ��Ϊһ����    �˶�ʱ����߳�����ԭ�߳��ܶ�       
%�о���������ȡ̫����˲��ʺϴ������� ����Ȼ�����ٶȴ�СҲ�й�ϵ

% Tension_calculation.m     ���������ڵ������������
% Tension_in_node.m         ����ǣ�����Խڵ������
% Tension_in_thread.m       �����������ߺ����߶ε�����

clear;
clc;
close all;
% kn1 kn2 �Ǽ���ǣ����������ϵ���� k1 k2�Ǽ��������߶ε�����ϵ�� 
global kn1 kn2 k1 k2 N dL0 mi mn

L = 5;                  % �����ܳ���
dL0 = 1;                % ������Ԫ����
Lct0 = 1.414*dL0;       % ǣ��������
N = L/dL0+1;            % �ڵ���  6��6
delta_t = 0.0001;       % ʱ�䲽��Ӱ��ܴ󣬲���̫�����ײ���̫��������0.0002�����Ӿͷ�ɢ�ˣ�0.0001����
tf = 1;                 % 40s       
Ns = round(tf/delta_t+1);      % �ܲ���
xnode = cell(Ns,1);     % cell����һ�㱻����Ԫ�����飬����ÿ����Ԫ���Դ��治ͬ���������ͣ���������ֵ���ַ�������Ԫ������ȣ�������ѧ����c������Ľṹ��
xnode4 = cell(Ns,1);

rho_thread = 1440;                         % kg/m3
Enet = 7e10;                               % ����ģ��  70GPA     
epsi= 0.106;                               % ����ȣ�damping ration of the net's thread
r_tether = 0.0001;                         % ���Ӱ뾶
r_ct = 0.002;                              % ǣ�����뾶
mi = pi*r_tether^2*dL0*rho_thread;         % ��������  
k1 = Enet*pi*r_tether^2/dL0;               % ����ϵ����k=N/m  Enet=N/m^2   EA/L0
k2 = 2*epsi*sqrt(mi*k1);                   % ����ϵ����cij=2*epsi(zeta)*sqrt(m*k)

mn = pi*r_ct^2*Lct0*rho_thread;            % ǣ�����ε�����
kn1 = Enet*pi*r_ct^2/Lct0;
kn2 = 2*epsi*sqrt(mn*kn1);

m=0.5;                                      % 4��ǣ��������������
m_sate = m + mn/2;                          % update��+ǣ������ε�����

T = zeros(3*N,N);
n0 = 0;                                     % ������ٶ� sqrt(398600/42164^3)
Ao = [zeros(3,3),eye(3)
    3*n0^2,0,0,0,2*n0,0
    0,0,0,-2*n0,0,0
    0,0,-n0^2,0,0,0];                    % cw���� 
Bo = [zeros(3,3);eye(3)];                % ���ӽڵ㷽��  ���ڵ���3������
Bo1 = [zeros(3,3);eye(3)];              
dx = @(t,xi,u)(Ao*xi+Bo*u);              % u  ���� ���ٶ�
dx1 = @(t,xi,u)(Ao*xi+Bo1*u);



%% initial position of net node
xnode{1} = zeros(6*N,N);           % ͨ����ά���Ͼ����ʾ��ά����ͦ��
xL0 = [10;0;L];                    % yzƽ��ĵڶ����� ����   ȷ���������Ͻ������ʼλ�ã�a��
x0 = cell(N,N);

for i = 1:N
    for j = 1:N
        xnode{1}((i-1)*6+1:(i-1)*6+3,j) = [xL0(1);xL0(2)+0.9*dL0*(j-1);xL0(3)-0.9*dL0*(i-1)];  %0.9��ʾδ���չ��    % ��ʼλ��ȷ�����ٶȶ���0(�ڹ������ϵ��Ҳ����) ��һ��ƽ�棬x��������ͬ
%       xsave{i,j}(1:3,1) = [xL0(1);xL0(2)+dL0*(j-1);xL0(3)-dL0*(i-1)];
    end
end

% ��ҷ��������ʼλ���ٶ�
xnode4{1}=zeros(6,4);

v=0;  
% ��ҷ������ ��ʼֵ   ��ʼ�ٶ� 1m/s
xnode4{1}(:,1)=[xL0(1)+Lct0;xL0(2)+dL0*(1-1);xL0(3)-dL0*(1-1); 1;-v;v];
xnode4{1}(:,2)=[xL0(1)+Lct0;xL0(2)+dL0*(N-1);xL0(3)-dL0*(1-1); 1;v;v];
xnode4{1}(:,3)=[xL0(1)+Lct0;xL0(2)+dL0*(N-1);xL0(3)-dL0*(N-1); 1;v;-v];
xnode4{1}(:,4)=[xL0(1)+Lct0;xL0(2)+dL0*(1-1);xL0(3)-dL0*(N-1); 1;-v;-v];
%% ����
for k = 1:Ns-1
    
    T4 = Tension_in_node(xnode{k},xnode4{k});  %4����������ҷ������֮�����
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
            u = T(3*(i-1)+1:3*i,j)/m_net(i,j);                  %u ���ٶ���
            [t,y] = ode45(@(t,xi)dx(t,xi,u),[(k-1)*delta_t,k*delta_t],xi);
            xi = y(end,:)';
            xnode{k+1}((i-1)*6+1:i*6,j) = xi;
%             xsave{i,j}(:,k+1) = xi;
        end
    end
   
 %% 4����  
 
 
 %xnode42(k,:)=xnode4(:,2)';    %��ҷ������1ʱ������ͼ  ���Էŵ�����ѭ����{}��������������Ҳ��������
 
 for i=1:4
        xi4 = xnode4{k}(:,i);
        u = -T4(:,i)/m_sate;
        [t,y] = ode45(@(t,xi)dx1(t,xi,u),[(k-1)*delta_t,k*delta_t],xi4);
        xi4 = y(end,:)';
        xnode4{k+1}(:,i) = xi4;       
 end 
        
    k
    
end
% ��ͼ
figure('color',[1 1 1]);        %ͼ�α�����Ϊ��ɫ��
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

