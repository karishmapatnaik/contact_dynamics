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
alpha0 = 0.05;

%% time specs
t = 0.0;
dt = 0.0025;
tend = 1;

%% initial conditions
v = [0 0 4]';
x = 0;
y = 1;
q = [x y theta]';
tfs = 0:dt:tend;
qfs = zeros(length(tfs),6);
lambdafs = zeros(length(tfs),7);
thetafs = zeros(length(tfs),1);

z0 = zeros(7,1);

Q = z0;

%% main loop
i = 1;

while t <= tend

    % set qhat_prev = 0;
    qhat_prev = 0*q;
        
    % Find n and alpha0 at this approx qhat
    n = [[0 1 l*cos(theta)/2]' [0 1 -l*cos(theta)/2]' ];
    D = [[1 0 l*sin(theta)/2]' [1 0 -l*sin(theta)/2]'  [-1 0 -l*sin(theta)/2]' [-1 0 l*sin(theta)/2]' ];
    f = [q(2)+(l/2)*sin(theta)-rho q(2)-(l/2)*sin(theta)-rho]';
%     if i == 155 %% at i == 155, its above contact, 
%     % % such that in the next step it will be violating constraint
%         break
%     end
    alpha0 = n'*q - f;

    %%%%%%%%%%%%%%%%%% fixed point iteration %%%%%%%%%%%%%%%%%

    % approximate qhat 
    qhat = q + dt*v;
    thetahat = qhat(3);
    
    while norm(qhat - qhat_prev)  > 0.005
    
        % Find n and alpha0 at this approx qhat
        n = [[0 1 l*cos(thetahat)/2]' [0 1 -l*cos(thetahat)/2]' ];
        D = [[1 0 l*sin(thetahat)/2]' [-1 0 -l*sin(thetahat)/2]'  [1 0 -l*sin(thetahat)/2]' [-1 0 l*sin(thetahat)/2]' ];
        f = [qhat(2)+(l/2)*sin(thetahat)-rho qhat(2)-(l/2)*sin(thetahat)-rho]';
        alpha0 = n'*qhat - f; 
    
        % Solve LCP for this new n and alpha0
        A = [D'*M^(-1)*D D'*M^(-1)*n e;
             n'*M^(-1)*D n'*M^(-1)*n [0;0];
             -e' mu 0 0];
        b = [D'*(v + dt*M^(-1)*[0 g 0]');
             (n'*q-alpha0)/dt + n'*(v+dt*M^(-1)*[0 g 0]');
             0];

        % using lemke 
        [Q,err] = lemke(A,b,z0);
        beta = Q(1:4);
        lambda = Q(7);
        cn = Q(5:6);
%         lambda = LCP(A,b);
        
        % update previous value
        qhat_prev = qhat;

        % next estimate
        vhat = v + M^(-1)*(dt*n*Q(5:6) + dt*D*Q(1:4) + [0 g*dt 0]');
        qhat = q + dt*v; 
    end

    % now time step
    v_ = v;
    v = v + M^(-1)*(n*Q(5:6) + D*Q(1:4) + [0 g*dt 0]');
    q = q + dt/2*(v_+v); 

    % store full state
    qfs(i,:) = [q' v'];
    theta = q(3);
    lambdafs(i,:) = Q';
    thetafs(i) = rad2deg(theta);

    % update time
    t = t + dt;
    i = i + 1;

end

figure
hold on
plot(tfs,qfs(:,end),'r')
plot(tfs,qfs(:,4),'b')
plot(tfs,qfs(:,5),'m')
grid minor

figure(12)
axis([-0.3 0.3 -0.1 y+0.2])
hold on
grid on
shg
lhs = [qfs(:,1)-((l/2)*cos(qfs(:,3))) qfs(:,2)-((l/2)*sin(qfs(:,3)))];
rhs = [qfs(:,1)+((l/2)*cos(qfs(:,3))) qfs(:,2)+((l/2)*sin(qfs(:,3)))];
j = i-1;
for i = 1:5:j
    figure(12)
    plot([lhs(i,1),rhs(i,1)],[lhs(i,2),rhs(i,2)],'-r')
    plot(qfs(i,1),qfs(i,2),'*b')
    pause(0.1)
end
