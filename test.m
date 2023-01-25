clc
clear all

% physical properties
m = 1;
k = 20;
b = 5;

% initial and desired conditions
xd = [1 0];
x = 0;
xdot = 1;
dt = 0.01;
tend = 10;
t = 0;
i = 1;

% variable for storing values
n = 0:dt:tend;
xfs = zeros(length(n),2);

% controller gains
kx = 20;
kv = 10;

while t < tend
    % code for MPC should go here and replace the controller
    
    % control input
    u = kx*(xd(1)-x) + kv*(xd(2)-xdot);
    % propagate states
    xdot = xdot + (-b*xdot+u)*dt;
    x = x + xdot*dt;
    % store for plotting
    xfs(i,:) = [x xdot];
    % increment time, counter
    t = t + dt;
    i = i + 1;
end


