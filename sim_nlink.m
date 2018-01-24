close all

m = 1;
L = 2;
g = 10;

tspan = 0:0.1:10;
y0 = [pi 0 0.01 0 0 0].'; %initial values for [theta1 theta2 theta3 omega1 omega2 omega3]
[t,y] = ode45(@(t,y)nlink(y,m,L,g,[0 0 0 0 0 0]),tspan,y0);%nlink(y,m,l,g,u -> input force vector

for k=1:length(t)
    drawmanip(y(k,:),m,L)
end