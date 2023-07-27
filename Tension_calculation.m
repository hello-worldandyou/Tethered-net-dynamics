%% 这程序编的贼棒，力的方向很准确
function [T,m_net] = Tension_calculation(i,j,T1,T2,T4)
% Given the index of the net node, calculating the tension and mass suffered by the node

global N mi mn

if i == 1 && j == 1
    node_type = 'a';
elseif i == 1 && j == N
    node_type = 'b';
elseif i == N && j == N
    node_type = 'c';
elseif i == N && j == 1
    node_type = 'd';                 %a b c d 绳网四个顶点
elseif i == 1 && j ~= 1 && j ~= N    
    node_type = 'r1';
elseif i == N && j ~= 1 && j ~= N
    node_type = 'rN';
elseif j == 1 && i ~= 1 && i ~= N
    node_type = 'c1';
elseif j == N && i ~= 1 && i ~= N
    node_type = 'cN';                %绳网除顶点外的四个边
else
    node_type = 'N';                 %内部点
end

testtension = mi*[0;0;0];
% testtension = mi*5*[0;0;0];        %作用在绳网顶点

switch node_type
    case 'a'    % The tension suffered by the node a
        T = T1(1:3,1)+T2(1:3,1)+testtension+T4(1:3,1);   %T4航天器与网角之间的力
        m_net = mi + mn/2;
    case 'b'
        T = -T1(1:3,end)+T2(1:3,end)+testtension+T4(1:3,2);
        m_net = mi + mn/2;
    case 'c'
        T = -T1(end-2:end,end)-T2(end-2:end,end)+testtension+T4(1:3,3);
        m_net = mi + mn/2;
    case 'd'
        T = T1(end-2:end,1)-T2(end-2:end,1)+testtension+T4(1:3,4);
        m_net = mi + mn/2;
    case 'r1'
        T = -T1(1:3,j-1)+T1(1:3,j)+T2(1:3,j);
        m_net = mi*1.5;
    case 'rN'
        T = -T1(end-2:end,j-1)+T1(end-2:end,j)-T2(end-2:end,j);
        m_net = mi*1.5;
    case 'c1'
        T = T1(3*(i-1)+1:3*i,1)-T2(3*(i-2)+1:3*(i-1),1)+T2(3*(i-1)+1:3*i,1);
        m_net = mi*1.5;
    case 'cN'
        T = -T1(3*(i-1)+1:3*i,end)-T2(3*(i-2)+1:3*(i-1),end)+T2(3*(i-1)+1:3*i,end);
        m_net = mi*1.5;
    case 'N'
        T = -T1(3*(i-1)+1:3*i,j-1)+T1(3*(i-1)+1:3*i,j)-T2(3*(i-2)+1:3*(i-1),j)+T2(3*(i-1)+1:3*i,j);
        m_net = mi*2;
end
