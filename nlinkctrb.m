syms m1 m2 m3 m l1 l2 l3 l1_cg l2_cg l3_cg th1 th2 th3 d t g w1 w2 w3 w1dot w2dot w3dot T1 T2 T3 sw

%m = 1; %1kg
% M = 5;
%l1 = 2; %2m

%%Complex Mass Distribution where successive mass is 1/2 the previous mass.
m1 = m;
m2 = m1/2;
m3 = m2/2;

%s = 1; % s = 1 pendulum is down, s = -1 pendulum is up

%%length distribution assuming constant density
l2 = l1/2;
l3 = l2/2;

%%every joint is at the center of oscillation
l1_cg = 2*l1/3;
l2_cg = 2*l2/3;
l3_cg = 2*l3/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th12 = th1 + th2;
th123 = th12 + th3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Coordinates of end points%%%
Pc_1 = [l1_cg*sin(th1); sw*l1_cg*cos(th1); 0];
Pc_2 = [l1_cg*sin(th1)+l2_cg*sin(th12); sw*(l1_cg*cos(th1)+l2_cg*cos(th12)); 0];
Pc_3 = [l1_cg*sin(th1)+l2_cg*sin(th12)+l3_cg*sin(th123); sw*(l1_cg*cos(th1)+l2_cg*cos(th12)+l3_cg*cos(th123)); 0];

% State Vector
q = [th1 th2 th3].';

%Linear Velocity Jacobians
Jv_1 = sym(zeros(3));
Jv_2 = sym(zeros(3));
Jv_3 = sym(zeros(3));

for i = 1:3
    for j = 1:3
        Jv_1(i,j) = diff(Pc_1(i,1), q(j,1));
        Jv_2(i,j) = diff(Pc_2(i,1), q(j,1));
        Jv_3(i,j) = diff(Pc_3(i,1), q(j,1));
    end
end

%Rotational Velocity Jacobians
Jw_1 = [0 0 0;
        0 0 0;
        1 0 0];     %all revolute joints
    
Jw_2 = [0 0 0;
        0 0 0; 
        1 1 0];
    
Jw_3 = [0 0 0;
        0 0 0;
        1 1 1];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inertia Matrices

Izz_1 = (m1*l1^2)/3;
Izz_2 = (m2*l2^2)/3;
Izz_3 = (m3*l3^2)/3;

Ic_1 = Izz_1;
Ic_2 = Izz_2;
Ic_3 = Izz_3;   %all revolute joints
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding mass matrix
%Due to linear velocities
m1_mat = m1*((Jv_1.')*Jv_1);
m2_mat = m2*((Jv_2.')*Jv_2);
m3_mat = m3*((Jv_3.')*Jv_3);
M_mat_Lin = simplify(m1_mat + m2_mat + m3_mat);

%Due to rotation velocities
I1_mat = ((Jw_1.')*Ic_1*Jw_1);
I2_mat = ((Jw_2.')*Ic_2*Jw_2);
I3_mat = ((Jw_3.')*Ic_3*Jw_3);
M_mat_Rot = I1_mat + I2_mat + I3_mat;

%Total Mass Matrix
M_mat = simplify(M_mat_Lin + M_mat_Rot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_i = [th1, th2, th3];
qqdot = [w1*w2 w1*w3 w2*w3].';
qdot_e2 = [w1^2 w2^2 w3^2].';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Christoffel Symbols
%b_ijk = (1/2)*(m_ikj + m_jki - m_ijk)

%Centrifugal C matrix
for i=1:3
    for j=1:3
        k=j;
        m_ijk = diff(M_mat(i,j),q_i(1,k));
        m_ikj = diff(M_mat(i,k),q_i(1,j));
        m_jki = diff(M_mat(j,k),q_i(1,i));
        
        Test1(i,j) = i*100+j*10+k;
        C_q(i,j) = simplify(0.5*(m_ikj + m_jki - m_ijk));
    end
end

%Coriolis B Matrix
count = 1;
k=0;
i=0;
j=0;
for i=1:3
    for kk=1:3
        k = kk+1;
        if k == 4
           k = 3;
        end
        %while k<=3
            if count<3
                j=k-count;
            else
                j=k-1;
                count = 0;
            end
            M_ijk = diff(M_mat(i,j),q_i(1,k));
            M_ikj = diff(M_mat(i,k),q_i(1,j));
            M_jki = diff(M_mat(j,k),q_i(1,i));
            B_q(i,kk) = simplify(M_ikj + M_jki - M_ijk);
            
            Test(i,kk) = i*100+j*10+k;
            
            count = count+1;
            k = k+1;
        %end
        
    end
    j=0;
end

%Gravity vector
g_v = [0 -g 0].';

%Gravity matrix (Potential Energy)
G_q = simplify(-Jv_1.'*m1*g_v-Jv_2.'*m2*g_v-Jv_3.'*m3*g_v);

V_q = B_q*qqdot + C_q*qdot_e2;

%Torque vector
Tau = [T1 T2 T3].';

V = subs(V_q, sw, 1);
G = subs(G_q, sw, 1);
M = subs(M_mat, sw, 1);

% Tau = sym(Mass_mat*qddot + V_q + G);

qddot = sym(simplify((M\(Tau - V - G)))); %qddot equations of motion

% Init Jacobian
J = sym(zeros(6));

% Combine all state equations into dx/dt
x_dot = simplify([w1 ;w2 ;w3 ;qddot]);

% Combine all state vectors into x (q)
q_fin = [q; w1; w2; w3];

% Method of finding Jacobian, differentiate each element of matrix with
% corresponding state vector
for i = 1:6
    for j = 1:6
        J(i,j) = diff(x_dot(i,1), q_fin(j,1));
    end
end

% Substitute 0 about all fixed points of stability
% Input forces tend to 0
% first order time derivatives tend to 0
% theta1, theta2, theta3 are all stable at intervals of 2n*pi()
A = subs(J, q_fin.', [0,0,0,0,0,0]);

BT = sym(M\Tau);

B1 = subs(BT, q_fin.', [0, 0, 0, 0, 0, 0]);

B = [0; 0; 0; B1];

% x_dot = A*x+ B*u

A
B