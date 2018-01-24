
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Can change values for different values of mass
%All successive link masses are in fractions of link1 mass
%refer to nlinkctrb for derivation
m = 1;
L = 2;
g = 10;

s = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_4_1 = s*(813*g)/(872*L);
A_4_2 = -s*(171*g)/(436*L);
A_4_3 = -s*(9*g)/(218*L);
A_5_1 = -s*(201*g)/(872*L);
A_5_2 = s*(1191*g)/(436*L);
A_5_3 = -s*(75*g)/(218*L);
A_6_1 = -s*(45*g)/(218*L);
A_6_2 = -s*(75*g)/(109*L);
A_6_3 = s*(474*g)/(109*L);

A = [
        0     0     0 1 0 0;
        0     0     0 0 1 0;
        0     0     0 0 0 1;
    A_4_1 A_4_2 A_4_3 0 0 0;
    A_5_1 A_5_2 A_5_3 0 0 0;
    A_6_1 A_6_2 A_6_3 0 0 0
    ];

B_4_1 = (9*(55)/(436*L^2*m));
B_4_2 = -(9*(123)/(436*L^2*m));
B_4_3 = (9*(20)/(436*L^2*m));
B_5_1 = -(9*(123))/(436*L^2*m);
B_5_2 = (9*(719))/(436*L^2*m);
B_5_3 = -(9*(996))/(436*L^2*m);
B_6_1 = (9*(5))/(109*L^2*m);
B_6_2 = -(9*(249))/(109*L^2*m);
B_6_3 = (9*(1508))/(109*L^2*m);

B = [
    0     0     0;
    0     0     0;
    0     0     0;
    B_4_1 B_4_2 B_4_3;
    B_5_1 B_5_2 B_5_3;
    B_6_1 B_6_2 B_6_3
    ];
%To check for controllability
if rank(A) == rank(ctrb(A,B))
    display('The system is controllable');
end

C = eye(6);
sys = ss(A,B,C,0*B);

tspan = 0:.001:50;
if(s==-1)
    y0 = [0 0 0 0 0 0].';
    [t,yNL] = ode45(@(t,y)nlink(y,m,L,g,[0 0 0 0 0 0]),tspan,y0);
elseif(s==1)
    y0 = [pi pi pi 0 0 0].';
    [t,yNL] = ode45(@(t,y)nlink(y,m,L,g,[0 0 0 0 0 0]),tspan,y0);
else    
end

figure
for k=1:100:length(t)
    drawmanip(yNL(k,:),m,L)
end