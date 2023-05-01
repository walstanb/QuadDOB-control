clear;
clc;

% skew-symmetric operator
ssmat = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
% inverse of skew-symmetric operator
vmap = @(x) [-x(2,3) x(1,3) -x(1,2)]';

% Vehicle properties
m = 4.34;
J = diag([0.0820, 0.0845, 0.1377]);

% Misc things
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';
g = 9.8;

% Initial conditions
x = [0 0 0]';
v = [0 0 0]';
R = eul2rotm([0 pi-0.01 0]);
omega = [0 0 0]';
Rd = zeros(3);
Rd_prev = zeros(3);

% Setpoint
xd = [0 0 10]';
vd = [0 0 0]';
b1d = [1 0 0]';

% Controller gains
kx = 10;
kv = 10;
kR = 10;
komega = 10;

% Main loop
dt = 0.001;
end_time = 10;

xlist = [];
ylist = [];
zlist = [];


for i = 1:end_time/dt
    % Translational error computations
    ex = x - xd;
    ev = v - vd;

    % Compute desired rotation
    b3d = -(-kx*ex - kv*ev - m*g*e3)/vecnorm(-kx*ex - kv*ev - m*g*e3);
    b2d = cross(b3d, b1d)/vecnorm(cross(b3d, b1d));
    Rd = [cross(b2d, b3d), b2d, b3d];

    % Compute desired angular velocity
    dRd = (Rd - Rd_prev)/dt;
    omegad = vmap(Rd'*dRd);
    Rd_prev = Rd;

    % Rotational/Angular error computations
    eR = 0.5*vmap(Rd'*R-R'*Rd);
    eomega = omega - R' * Rd * omegad;

    % Control input computation
    f = -dot((-kx*ex - kv*ev - m*g*e3), R*e3);
    M = -kR*eR - komega*eomega + cross(omega, J*omega) - J*(ssmat(omega)*R'*Rd*omegad);

    % Control Allocation Matrix
    u = [f M']';
    d = 0.25; ctf = 0.1;
    CAM = [1 1 1 1; 0 -d 0 d; d 0 -d 0; -ctf ctf -ctf ctf];
    fi = CAM\u;
    fi = min(fi, 20);
    fi = max(fi, -20);
    u = CAM * fi;
    f = u(1);
    M = u(2:4);

    % System evolution
    a = g*e3 - f*R*e3/m;
    v = v + a*dt;
    x = x + v*dt;
    
    xlist(i) = x(1);
    ylist(i) = x(2);
    zlist(i) = x(3);

    domega = J\(M - cross(omega,J*omega));
    omega = omega + domega * dt;
    R = R + R * ssmat(omega) * dt;
    R = orthnorm(R);
    
    aHistory(:,i) = a;
    vHistory(:,i) = v;
    xHistory(:,i) = x;
    domegaHistory(:,i) = domega;
    omegaHistory(:,i) = omega;
    RHistory(:,:,i) = R;
    exHistory(:,i) = ex;
    evHistory(:,i) = ev;
    eRHistory(:,i) = eR;
    eomegaHistory(:,i) = eomega;
end

figure
plot3(xlist, ylist, zlist);



% time = 0:dt:end_time-dt;
% figure(1);clf
% sp331 = subplot(3,3,1); hold on; grid on; grid minor;
% plot(time, aHistory); legend("ax", "ay", "az"); title("Linear Acceleration");
% sp332 = subplot(3,3,2); hold on; grid on; grid minor;
% plot(time, vHistory); legend("vx", "vy", "vz"); title("Linear Velocity");
% sp333 = subplot(3,3,3); hold on; grid on; grid minor;
% plot(time, xHistory); legend("x", "y", "z"); title("Linear Position");
% sp334 = subplot(3,3,4); hold on; grid on; grid minor;
% plot(time, domegaHistory); legend("dp", "dq", "dr"); title("Angular Acceleration");
% sp335 = subplot(3,3,5); hold on; grid on; grid minor;
% plot(time, omegaHistory); legend("p", "q", "r"); title("Angular Velocity");
% sp336 = subplot(3,3,6); hold on; grid on; grid minor;
% plot(time, rotm2eul(RHistory)); legend("yaw", "pitch", "roll"); title("Euler Angles");
% sp337 = subplot(3,3,7); hold on; grid on; grid minor;
% plot(time, exHistory); legend("ex", "ey", "ez"); title("Position Error");
% sp338 = subplot(3,3,8); hold on; grid on; grid minor;
% plot(time, eRHistory); legend("eRx", "eRy", "eRz"); title("Rotational Error");
% sp339 = subplot(3,3,9); hold on; grid on; grid minor;
% plot(time, eomegaHistory); legend("eomegax", "eomegay", "eomegaz"); title("Anguler Velocity Error");
% 

function R = orthnorm(R)
    R(:,1) = R(:,1)/vecnorm(R(:,1));
    R(:,2) = cross(R(:,3), R(:,1))/vecnorm(cross(R(:,3), R(:,1)));
    R(:,3) = cross(R(:,1), R(:,2))/vecnorm(cross(R(:,1), R(:,2)));
end