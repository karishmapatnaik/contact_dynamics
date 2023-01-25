clc
clear

% physical properties
% k = 850;
% b = 6.5;
% m = 0.0015;
m = 0.0015;
k = 850*m;
b = 6.5*m;
A = [0 1 0; 0 0 0; 0 0 -k/b];
B = [0 ; 1/m; -1/b];
C = [1 0 -1];
z = 3; % position of mass
l = 0.1; % distance between spring end and cog when spring at rest
x = z - l; % change in l when compressed 
zdot = -20; % initial condition when spring at rest
xdot = zdot;
zeta = l; % when spring at rest, this is same as l
zeta_bar = zeta - l;
eta = [x; xdot; zeta_bar];
etadot = [xdot; 0; 0];

% simulation specs
tend = 0.1;
t = 0;
dt = 0.0025;
tfs = (0:dt:tend)';
etafs = zeros(length(tfs),6);
etafs(1,:) = [eta' etadot'];
zfs = zeros(length(tfs),1);
lambdafs = zeros(length(tfs),1);
i = 1;

% Time stepping
while t < tend+1
    % LCP
    [lambda,err] = lemke(C*B,C*A*eta,0);
    % time-stepping
    etadot = A*eta + B*lambda;
    eta = eta + etadot*dt;
    % storing full state
    etafs(i+1,:) = [eta' etadot'];
    lambdafs(i) = lambda;
    zfs(i) = eta(1) + l;
    % update time
    t = t + dt;
    i = i + 1;
end

plot(tfs,zfs(1:length(tfs),1))
shg
grid minor


