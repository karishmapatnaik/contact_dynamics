clc
clear

dt = 0.008; %% 0.001
tend = 100;
k = 8;
b = 0.1;

%% validation using step input
J = 1/(650); %% tuned to get an idea of inertia

%  diffferent initial conditions
tau = 0; %% for 6N force the angle should deflect 33^o
vspring = 1.5;

% % using force
% tau = 6*(113/1000); %% for 6N force the angle should deflect 33^o
% vspring = 0;

xspring = 0;
t = 0;
xfs_v = [xspring vspring];

for i = 1:tend/dt
    vspring = vspring + (tau/J -b*vspring - k*xspring)*dt;
    xspring = xspring + vspring*dt;
    xfs_v = [xfs_v; xspring vspring];
    t = t + dt;
end

tfs = linspace(0,t,length(xfs_v));

figure(10)
hold on
plot(tfs,xfs_v(:,1),'LineWidth',2)

rad2deg(xspring);

b^2 - 4*k*J