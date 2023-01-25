clc
clear

% physical properties
m = 0.0015;
k = 80*m;
b = 0.1*m;
z = 1; % position of mass
l = 0.08; % distance between spring end and cog when spring at rest
x = z - l; % change in l when compressed 
zdot = -0.1; % initial condition when spring at rest
xdot = zdot;
zeta = l; % when spring at rest, this is same as l

% simulation specs
tend = 16;
t = 0;
dt = 0.0025;
tfs = (0:dt:tend)';
zfs = zeros(length(tfs),3);
lambdafs = zeros(length(tfs),1);
zfs(1,:) = [z x xdot];
i = 2;

% Time stepping
while t < tend
    % LCP
    [lambda,err] = lemke(1,k*x+b*xdot,0);
    % time-stepping
    xdot = xdot + lambda/m*dt;
    x = x + xdot*dt;
    % storing full state = ;
    z = x + l;
    zfs(i,:) = [z x xdot];
    lambdafs(i) = lambda;
    % time updation
    t = t + dt;
    i = i + 1;
end

plot(tfs,zfs(1:length(tfs),1))
hold on
% plot(tfs,zfs(1:length(tfs),2))
grid on
grid minor

