clc
clear
close all

%% rigid body specifications
g = -9.81;
m = 1;
J = 0.002;
mu = 0.6;
M = [m 0 0; 0 m 0; 0 0 J];
e = [1 1 1 1]';
theta = pi/6;
l = 0.5;
rho = 0.05;

%% time specs
t = 0.0;
dt = 0.0025;
tend = 4;

%% initial conditions
v = [0 0 4]';
x = 0;
y = 0.5;
q = [x y theta]';
tfs = 0:dt:tend;
qfs = zeros(length(tfs),6);
lambdafs = zeros(length(tfs),7);

z0 = zeros(7,1);

%% main loop
i = 1;
while t <= tend
    n = [[0 1 l*cos(theta)/2]' [0 1 -l*cos(theta)/2]' ];
    D = [[1 0 l*sin(theta)/2]' [-1 0 l*sin(theta)/2]'  [-1 0 -l*sin(theta)/2]' [1 0 -l*sin(theta)/2]' ];
    f = [q(2)+(l/2)*sin(theta)-rho q(2)-(l/2)*sin(theta)-rho]';
    alpha0 = n'*q - f; 
%    solve LCP
    A = [D'*M^(-1)*D D'*M^(-1)*n e;
         n'*M^(-1)*D n'*M^(-1)*n [0;0];
         -e' mu 0 0];
    b = [D'*(v + dt*M^(-1)*[0 g 0]');
         (n'*q-alpha0)/dt + n'*(v+dt*M^(-1)*[0 g 0]');
         0];
    % using lemke
    [lambda,err] = lemke(A,b,z0);
% %     % using gauss-newton
%     lambda = LCP(A,b,[],[],z0,1);
    % time-stepping dynamics
    v_ = v;
    v = v + M^(-1)*(n*lambda(5:6) + D*lambda(1:4) + [0 g*dt 0]');
    q = q + dt/2*(v + v_);
    % store full state
    qfs(i,:) = [q' v'];
    theta = q(3);
    lambdafs(i,:) = lambda';
    % update time
    t = t +dt;
    i = i +1;
%     if i == 150 %% debugging what happens at i == 76, why is lambda(7) = 0.6
%         break
%     end
end

figure(10)
plot(tfs,qfs(:,4),'r')
axis([0 tend -0.5 0.5])
ylabel('Horizontal Velocity')
grid minor

figure(11)
plot(tfs,qfs(:,5),'b')
axis([0 tend -3 2])
ylabel('Vertical Velocity')
grid minor

figure(13)
plot(tfs,qfs(:,6),'m')
axis([0 tend -20 21])
ylabel('Angular Velocity')
xlabel('Time in (s)')
grid minor

figure(100)
axis([-0.3 0.3 -0.1 y+0.2])
hold on
grid on
shg
lhs = [qfs(:,1)-((l/2)*cos(qfs(:,3))) qfs(:,2)-((l/2)*sin(qfs(:,3)))];
rhs = [qfs(:,1)+((l/2)*cos(qfs(:,3))) qfs(:,2)+((l/2)*sin(qfs(:,3)))];
j = i-1;
for i = 1:1:j
    figure(100)
    plot([lhs(i,1),rhs(i,1)],[lhs(i,2),rhs(i,2)],'-r')
    N(i) = getframe;
    pause(0.1)
end
