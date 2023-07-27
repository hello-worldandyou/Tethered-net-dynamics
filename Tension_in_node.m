function T4 = Tension_in_node(xn,xn4)
% xn: �������нڵ��״̬
% xn4: 4��ǣ���ڵ��״̬
% N: ����ÿһ�ߵĽڵ���
% dL0: ����С����߳�

global N dL0 kn1 kn2

T4 = zeros(3,4); % the tension of the four drag thread
x1=zeros(6,4);
x1(:,1) = xn(1:6,1);
x1(:,2) = xn(1:6,N);
x1(:,3) = xn(6*(N-1)+1:6*N,N);
x1(:,4) = xn(6*(N-1)+1:6*N,1);
x2 = xn4;

dL04=1.414*dL0;
%% Tension in thread
% Tension of the thread in the row of the net
for i = 1:4
        
        L_t = -x1(1:3,i)+x2(1:3,i);         %�����ǽڵ�ָ������
        e_Lt = L_t/norm(L_t);               %��λ�� �����������λ��
        Lt = norm(L_t);                     
        Delta_v = -x1(4:6)+x2(4:6);         %�����ǽڵ�ָ������
        T4_norm = kn1*(Lt-dL04)+kn2*dot(Delta_v,e_Lt);
        if Lt < dL04
            T4(:,i) = zeros(3,1);
% �޸ĺ�
%         elseif Lt == dL0 && Delta_v'*e_Lt > 0  
%             T1(3*(i-1)+1:3*i,j) = k2*(Delta_v'*e_Lt)*e_Lt;
        elseif Lt > dL04
        %elseif (Lt > dL04) && (T4_norm>0)
            T4(:,i) = T4_norm * e_Lt;   %ע�����ķ���  չ������Ӧ�����и���С������
        else
            T4(:,i) = zeros(3,1);
        end
    end
