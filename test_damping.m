clc
clear 
close all

m = 0.0015;
b = 6.5*m;

dt = 0.0025;
tend = 10;
t = 0;

xdot = -20;
x = 3;
i = 1;

while t <= tend
    xdot = xdot + (-b*xdot)*dt/m;
    x = x + xdot*dt;
    t = t + dt;
    xfs(i,:) = [x xdot];
    i = i + 1;
end

plot(linspace(0,tend,length(xfs)),xfs(:,1))