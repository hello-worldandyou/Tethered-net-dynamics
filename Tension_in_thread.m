function [T1,T2] = Tension_in_thread(x)
% Given the orbit state vector of the net nodes in the orbit reference
% frame, calculating the tension in the thread of the net

% x: the orbit state vector of the net node
% T1: the tension of the thread in the row of the net
% T2: the tension of the thread in the column of the net
% N: 绳网每一边的节点数
% dL0: 绳网小网格边长
global k1 k2 N dL0

T1 = zeros(3*(N),N-1); % the tension of the thread in the row of the net  画个九宫格就明白了
T2 = zeros(3*(N-1),N); % the tension of the thread in the column of the net

x1 = x;
x2 = x;

%% Tension in thread
% Tension of the thread in the row of the net
for i = 1:N
    for j = 1:N-1
        xx1 = x1((i-1)*6+1:i*6,j);
        xx2 = x1((i-1)*6+1:i*6,j+1);
        L_t = -xx1(1:3)+xx2(1:3);           %同行相邻两列相对位置   向量→
        Lt = norm(L_t); 
        e_Lt = L_t/Lt;               %单位化 相邻两列相对位置
                                            %同行相邻点的距离
        Delta_v = -xx1(4:6)+xx2(4:6);       %相邻点相对速度        速度→
        T1_norm = k1*(Lt-dL0)+k2*dot(Delta_v,e_Lt);
        if Lt < dL0
            T1(3*(i-1)+1:3*i,j) = zeros(3,1);  
        elseif (Lt > dL0) && (T1_norm>0)
            T1(3*(i-1)+1:3*i,j) = T1_norm *e_Lt;
        else 
            T1(3*(i-1)+1:3*i,j) = zeros(3,1);                           
        end
    end
end

% Tension of the thread in the column of the net
for i = 1:N-1
    for j = 1:N
        xx1 = x2((i-1)*6+1:i*6,j);
        xx2 = x2(i*6+1:(i+1)*6,j);       % 方向 ↓
        L_t = -xx1(1:3)+xx2(1:3);
        Lt = norm(L_t);
        e_Lt = L_t/Lt;
        T2_norm = k1*(Lt-dL0)+k2*dot(Delta_v,e_Lt);
        Delta_v = -xx1(4:6)+xx2(4:6);    % 方向 ↓
        if Lt < dL0
            T2(3*(i-1)+1:3*i,j) = zeros(3,1);
        %elseif Lt == dL0 && Delta_v'*e_Lt > 0
         %   T2(3*(i-1)+1:3*i,j) =1*e_Lt; % k2*(Delta_v'*e_Lt)*e_Lt;
        elseif (Lt > dL0)
        %elseif (Lt > dL0) && (T2_norm>0)
             T1(3*(i-1)+1:3*i,j) = T2_norm *e_Lt;
        else 
            T2(3*(i-1)+1:3*i,j) = zeros(3,1);
        end
    end
end