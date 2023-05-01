clear

% Declaring Globals
global e3 R W  x v xt yt zt xdt ydt zdt xddt yddt zddt dt time;
global dbarx dbarW timesteps Rdes_previous wdes_previous I Rdes;
global zW zx lW lv Km Ks fb;
animate=0;

%Initialising all Variables
J = diag([0.0820, 0.0845, 0.1377]);
d = 0.315;
ctf = 8.004e-4;
m = 4.34;
g = 9.81;
upper_limit = 60;

kx = 16*m;
kv = 5.6*m;
kR = 8.81;
kw = 2.54;
Km = diag([0.2 0.2 0.05]);
Ks = 0.02;

dt = 0.01;
time = 20;

reset_params();

desired_path();
global xdest b1dest vdest adest;

for n = 1 : length(timesteps)

    % desired values
    xdes = xdest(:,n);
    b1des = b1dest(:,n); 
    vdes = vdest(:,n);
    ades = adest(:,n);

    % errors
    xerr = x - xdes;
    verr = v - vdes;
    
    xerr_hist(n) = xerr(1,1);
    yerr_hist(n) = xerr(2,1);
    zerr_hist(n) = xerr(3,1);

    b3des = (-kx*xerr - kv*verr - m*g*e3 + m*ades)/norm(-kx*xerr - kv*verr - m*g*e3 + m*ades);
    b2des = (cross(b3des, b1des)/norm(cross(b3des, b1des)));
    
    Rdes = [cross(b2des, b3des) b2des b3des];

    dRdes = (Rdes - Rdes_previous)/dt;
    Rdes_previous = Rdes;
    wdes = Vhelper(Rdes'*dRdes);

    werr = W - R'*Rdes*wdes;
    rerr = Vhelper(0.5*(Rdes'*R - R'*Rdes));

    % control inputs
    f = -dot((-kx*xerr - kv*verr-m*g*e3 + m*ades), R*e3);
    f = min(max(abs(f)), upper_limit)* sign(f);
    M = -kR*rerr - kw*werr;

%     f = min(max(abs(f)), upper_limit)/max(abs(f))*f;
%     M = min(max(abs(M)), 1)/max(abs(M))*M

    % Control Allocation Matrix
%     u = [f M']';
% 
%     CAM = [1 1 1 1; 0 -d 0 d; d 0 -d 0; -ctf ctf -ctf ctf];
%     fi = CAM\u;
%     
%     u = CAM * fi;
%     f = u(1)
%     M = u(2:4);

    vdot = (m*g*e3 - f*R*e3)/m;
    v = v + vdot*dt;
%     v = min(max(abs(v)), 1.4)/max(abs(v))*v;
    xdot = v;
    x = x + xdot*dt;
    Wdot = inv(J) * (M - cross(W,J * W));
    W = W + Wdot*dt;
    Rdot = R * Whelper(W);
    R = R + Rdot*dt;
    R = GSnorm(R);

    xt(n) = x(1);
    yt(n) = x(2);
    zt(n) = x(3);
   
    xdt(n) = v(1, 1);
    ydt(n) = v(2, 1);
    zdt(n) = v(3, 1);
    
    xddt(n) = (xdt(end)-xdt(end-1))/dt;
    yddt(n) = (ydt(end)-ydt(end-1))/dt;
    zddt(n) = (zdt(end)-zdt(end-1))/dt;
    
end


ax = axes(figure(1));
plot3(ax, xt,yt,zt);
xlabel('X');
ylabel('Y');
zlabel('Z');
hold(ax, "on");

axs = [];
val = [xt; yt; zt; xdt; ydt; zdt; xddt; yddt; zddt];
ylabels = ["x" "y" "z" "v_x" "v_y" "v_z" "a_x" "a_y" "a_z"];
fig2 = figure(2);
for i = 1 : 9
    axs(i) = subplot(3, 3, i);
    plot(axs(i), timesteps, val(i,:))
    xlabel('time') 
    ylabel(ylabels(i))
    hold(axs(i), "on");
end


erraxs = [];
val2 = [xerr_hist; yerr_hist; zerr_hist];
ylabels2 = ["Error in x" "Error in y" "Error in z"];
fig3 = figure(3);
for i = 1 : 3
    erraxs(i) = subplot(1, 3, i);
    plot(erraxs(i), timesteps, val2(i,:))
    xlabel('time') 
    ylabel(ylabels2(i))
    hold(erraxs(i), "on");
end

if animate
    scenario = uavScenario("StopTime", 8, "UpdateRate", 1000);
    axnm = axes(figure("Units","normalized", "OuterPosition", [0 0 1 1]));
    hold(axnm, "on")
    
    plot3(axnm, 70*yt', 70*xt', -70*zt');
end


reset_params();

for n = 1 : length(timesteps)
    [dis_x, dis_W] = dis_gen(timesteps(n));

    % desired values
    xdes = xdest(:,n);
    b1des = b1dest(:,n); 
    vdes = vdest(:,n);
    ades = adest(:,n);

    % errors
    xerr = x - xdes;
    verr = v - vdes;

    xerr_hist(n) = xerr(1,1);
    yerr_hist(n) = xerr(2,1);
    zerr_hist(n) = xerr(3,1);

    b3des = (-kx*xerr - kv*verr - m*g*e3 + m*ades)/norm(-kx*xerr - kv*verr - m*g*e3 + m*ades);
    b2des = (cross(b3des, b1des)/norm(cross(b3des, b1des)));
    Rdes = [cross(b2des, b3des) b2des b3des];

    dRdes = (Rdes - Rdes_previous)/dt;
    Rdes_previous = Rdes;
    wdes = Vhelper(Rdes'*dRdes);

    werr = W - R'*Rdes*wdes;
    rerr = Vhelper(0.5*(Rdes'*R - R'*Rdes));

    % control inputs
    f = -dot((-kx*xerr - kv*verr-m*g*e3 + m*ades), R*e3);
    f = min(max(abs(f)), upper_limit)* sign(f);
    M = -kR*rerr - kw*werr;

    % Control Allocation Matrix
%     u = [f M']';
%     CAM = [1 1 1 1; 0 -d 0 d; d 0 -d 0; -ctf ctf -ctf ctf];
%     fi = CAM\u;
%     fi = min(fi, upper_limit);
%     fi = max(fi, lower_limit);
%     u = CAM * fi;
%     f = u(1);
%     M = u(2:4);

    vdot = (m*g*e3 - f*R*e3 + dis_x)/m;
%     vdot = min(max(abs(vdot)), 1.4)/max(abs(vdot))*vdot;
    v = v + vdot*dt;
    xdot = v;
    x = x + xdot*dt;
    Wdot = inv(J) * (M - cross(W,J * W)+ dis_W);
    W = W + Wdot*dt;
    Rdot = R * Whelper(W);
    R = R + Rdot*dt;
    R = GSnorm(R);

    xt(n) = x(1);
    yt(n) = x(2);
    zt(n) = x(3);
   
    xdt(n) = v(1, 1);
    ydt(n) = v(2, 1);
    zdt(n) = v(3, 1);
    
    xddt(n) = (xdt(end)-xdt(end-1))/dt;
    yddt(n) = (ydt(end)-ydt(end-1))/dt;
    zddt(n) = (zdt(end)-zdt(end-1))/dt;
end

plot3(ax, xt, yt, zt);

if animate
    plot3(axnm, 70*yt', 70*xt', -70*zt');
end

val = [xt; yt; zt; xdt; ydt; zdt; xddt; yddt; zddt];
for i = 1 : 9
    plot(axs(i), timesteps, val(i,:))
end

val2 = [xerr_hist; yerr_hist; zerr_hist];
for i = 1:3
    plot(erraxs(i), timesteps, val2(i,:))
end

reset_params();

fblog_x = [];
fblog_y = [];
fblog_z = [];
dbarxlog_x = [];
dbarxlog_y = [];
dbarxlog_z = [];

for n = 1 : length(timesteps)
    
    [dis_x, dis_W] = dis_gen(timesteps(n));

    % desired values
    xdes = xdest(:,n);
    b1des = b1dest(:,n); 
    vdes = vdest(:,n);
    ades = adest(:,n);

    xerr = x - xdes;
    verr = v - vdes;

    xerr_hist(n) = xerr(1,1);
    yerr_hist(n) = xerr(2,1);
    zerr_hist(n) = xerr(3,1);

    b3des = (-kx*xerr - kv*verr - m*g*e3 + m*ades)/norm(-kx*xerr - kv*verr - m*g*e3 + m*ades);
    b2des = (cross(b3des, b1des)/norm(cross(b3des, b1des)));
    
    Rdes = [cross(b2des, b3des) b2des b3des];

    dRdes = (Rdes - Rdes_previous)/dt;
    Rdes_previous = Rdes;
    wdes = Vhelper(Rdes'*dRdes);
    dwdes = (wdes - wdes_previous)/dt;
    wdes_previous = wdes; 

    werr = W - R'*Rdes*wdes;
    rerr = Vhelper(0.5*(Rdes'*R - R'*Rdes));
    
    vdot = (m*g*e3 - R*fb + dis_x)/m;
    v = v + vdot*dt;
    xdot = v;
%     v = min(max(abs(v)), 1.4)/max(abs(v))*v;
    x = x + xdot*dt;
    Wdot = inv(J) * (M - cross(W,J * W)+ dis_W);
    W = W + Wdot*dt;
    Rdot = R * Whelper(W);
    R = R + Rdot*dt;
    R = GSnorm(R);

    lamdaW = Km * wdes;
    lamdav = Ks * vdes;
    
    M = -kR*rerr - kw*werr + J*dwdes + cross(W, J*W) - dbarW;

    f_lim = -dot((-kx*xerr - kv*verr - m*g*e3 + m*ades - dbarx), R*e3);
    f_lim = min(max(abs(f_lim)), upper_limit)* sign(f_lim);
    fb = f_lim * e3;


    
    
    fblog_x(n) = fb(1,1);
    fblog_y(n) = fb(2,1);
    fblog_z(n) = fb(3,1);
    dbarxlog_x(n) = dbarx(1,1);
    dbarxlog_y(n) = dbarx(2,1);
    dbarxlog_z(n) = dbarx(3,1);

    zWdot = - lW * [inv(J)*(lamdaW + zW) - inv(J)*(cross(W, J*W)) + inv(J)*M];
    dbarW = zW + lamdaW;
    zW = zW + zWdot*dt;

    zxdot = -lv * (1/m * (lamdav + zx) - g*e3 + ((R*fb)/m));
    dbarx = zx + lamdav;
    zx = zx + zxdot*dt;
    
    xt(n) = x(1);
    yt(n) = x(2);
    zt(n) = x(3);
   
    xdt(n) = v(1, 1);
    ydt(n) = v(2, 1);
    zdt(n) = v(3, 1);
    
    xddt(n) = (xdt(end)-xdt(end-1))/dt;
    yddt(n) = (ydt(end)-ydt(end-1))/dt;
    zddt(n) = (zdt(end)-zdt(end-1))/dt;
    
end


plot3(ax, xt, yt, zt, '--k');

xdes_path = xdest(1, :);
y_desired_path = xdest(2, :);
z_desiredpath = xdest(3, :);
plot3(ax, xdes_path, y_desired_path, z_desiredpath)
legend(ax, "Drone in Ideal Conditions", "With Disturbance", "With DOB", "Desired")
% hold(ax, "off");

val = [xt; yt; zt; xdt; ydt; zdt; xddt; yddt; zddt];
for i = 1 : 9
    plot(axs(i), timesteps, val(i,:), '--k')
    hold(axs(i), "off");
end
Lgnd = legend("Drone in Ideal Conditons Ideal", "With Disturbance", "With DOB");
Lgnd.Position(1) = 0.84;
Lgnd.Position(2) = 0.88;

val2 = [xerr_hist; yerr_hist; zerr_hist];
for i = 1:3
    plot(erraxs(i), timesteps, val2(i,:), '--k')
    hold(erraxs(i), "off");
end

Lgnd = legend("Drone in Ideal Conditons Ideal", "With Disturbance", "With DOB");
Lgnd.Position(1) = 0.84;
Lgnd.Position(2) = 0.88;

% figure(4)
% 
% subplot(2,3,1);
% plot(timesteps, dbarxlog_x, '--k');
% 
% subplot(2,3,4);
% plot(timesteps, fblog_x);
% 
% subplot(2,3,2);
% plot(timesteps, dbarxlog_y, '--k');
% 
% subplot(2,3,5);
% plot(timesteps, fblog_y);
% 
% subplot(2,3,3);
% plot(timesteps, dbarxlog_z, '--k');
% 
% subplot(2,3,6);
% plot(timesteps, fblog_z);
% 
% legend("dbarx", "fb")

if animate
    scale = 1
    plot3(axnm, scale*yt', scale*xt', -scale*zt', '--k');
    
    plat = uavPlatform("UAV", scenario, "Trajectory", ...
        waypointTrajectory(70*[xt' yt' zt'], "TimeOfArrival", timesteps));
    updateMesh(plat,"quadrotor", {10}, [1 0 0], eul2tform([0 0 pi]));
    
    [axnm, plotFrames] = show3D(scenario);
    view([-110 30]);
    axis equal
    
    setup(scenario);
    updateCounter = 0;
    while true
        % Advance scenario.
        isRunning = advance(scenario);
        updateCounter = updateCounter + 1;
    
        % Update visualization every 10 updates.
        if updateCounter > 10
            show3D(scenario, "FastUpdate", true, "Parent", axnm);
            updateCounter = 0;
            drawnow limitrate
        end
    
        if ~isRunning 
            break; 
        end
    end
    hold(axnm, "off")
end

function Wh = Whelper(W)
    Wh = [0    -W(3)  W(2);
          W(3)  0    -W(1);
         -W(2)  W(1)   0];
end

function Vh = Vhelper(W)
    Vh = [W(3, 2); W(1, 3) ; W(2,1)];
end

function b = GSnorm(a) 
    [m,n]=size(a);
    b=zeros(m,n);
    k=n;
    b(:,1)=a(:,1)./sqrt(sum(a(:,1).*a(:,1))); 
    b(:,2)=a(:,2)-dot(b(:,1),a(:,2))*b(:,1);
    b(:,2)=b(:,2)./sqrt(sum(b(:,2).*b(:,2)));
    if k>2
        for i=3:k
            b(:,i)=a(:,i)-dot(b(:,2),a(:,i))*b(:,2)-dot(b(:,1),a(:,i))*b(:,1);
            if i>3
                for j=i-1:-1:3
                    b(:,i)=b(:,i)-dot(b(:,j),a(:,i))*b(:,j);
                end
                
            end
         b(:,i)=b(:,i)./sqrt(sum(b(:,i).*b(:,i)));
        end
    end
end

function reset_params()
    global e1 e2 e3 R W Wh x v xt yt zt xdt ydt zdt xddt yddt zddt dt time;
    global dbarx dbarW timesteps Rt Rdes_previous wdes_previous I Rdes;
    global zW zx lW lv Km Ks fb;

    e1 = [1; 0; 0];
    e2 = [0; 1; 0];
    e3 = [0; 0; 1];
    R = [e1 e2 e3];
    W = [0; 0; 0];
    Wh = [0    -W(3)  W(2);
          W(3)  0    -W(1);
         -W(2)  W(1)   0];
    x = [0 0 0]';
    v = [0 0 0]';
    
    xt = [];
    yt = [];
    zt = [];
    
    xdt = [0 0];
    ydt = [0 0];
    zdt = [0 0];
    
    xddt = [];
    yddt = [];
    zddt = [];
    
    dbarx = [0;0;0];
    dbarW = 0;
    
    xerr_hist = [];
    yerr_hist = [];
    zerr_hist = [];

    %The timestep calculations
    timesteps = 0:dt:time;
    
    Rt = cell(length(timesteps), 1);
    Rdes_previous = zeros(3);
    wdes_previous = [0; 0; 0];
    Rdes = zeros(3);
    
    I = eye(3);
    
    zW = [0; 0; 0];
    zx = [0; 0; 0];
    
    lW = Km * I;
    lv = Ks * I;
    
    fb = [0 0 0]';
end

function desired_path()

    global xdest b1dest vdest adest timesteps;
    
    xdest = [];
    b1dest = [];
    vdest = [];
    adest = [];

%     syms t;
%     xdfn = [0.4*t; 0.4*sin(pi*t); 0.6*cos(pi*t)];
% 
%     for n = 1 : length(timesteps)
%         xdest(:,n) = [0.4*timesteps(n) 0.4*sin(pi*timesteps(n)) 0.6*cos(pi*timesteps(n))]';
%         b1dest(:,n) = [cos(pi*timesteps(n)) sin(pi*timesteps(n)) 0]';
%         vdest(:,n) = [0.4 0.4*pi*cos(pi*timesteps(n)) -0.6*pi*sin(pi*timesteps(n))]';
%         adest(:,n) = [0 -0.4*pi^2*sin(pi*timesteps(n)) -0.6*pi^2*cos(pi*timesteps(n))]';
%     end
    
    for n = 1 : length(timesteps)
        xdest(:,n) = [0.6*cos(pi*timesteps(n)) 0.4*sin(pi*timesteps(n)) 0.6*timesteps(n)]';
        b1dest(:,n) = [0 sin(pi*timesteps(n)) cos(pi*timesteps(n))]';
        vdest(:,n) = [-0.6*pi*sin(pi*timesteps(n)) 0.4*pi*cos(pi*timesteps(n)) 0.6]';
        adest(:,n) = [-0.6*pi^2*cos(pi*timesteps(n)) -0.4*pi^2*sin(pi*timesteps(n)) 0]';
    end

%     syms t;
%     xdfn = [2*cos((pi*t)/10); 2*sin((pi*t)/10); (t/2)]
%     b1fn = diff(xdfn) / norm(diff(xdfn))
%     diff(xdfn)
%     diff(diff(xdfn))

%     for n = 1 : length(timesteps)
%         xdest(:,n) = [2*cos((pi*timesteps(n))/10); 2*sin((pi*timesteps(n))/10); (timesteps(n)/2)];
%         b1dest(:,n) = [-(pi*sin((pi*timesteps(n))/10))/(5*((pi^2*abs(cos((timesteps(n)*pi)/10))^2)/25 + (pi^2*abs(sin((timesteps(n)*pi)/10))^2)/25 + 1/4)^(1/2)); 
%             (pi*cos((pi*timesteps(n))/10))/(5((pi^2*abs(cos((timesteps(n)pi)/10))^2)/25 + (pi^2*abs(sin((timesteps(n)*pi)/10))^2)/25 + 1/4)^(1/2)); 
%             1/(2((pi^2*abs(cos((timesteps(n)*pi)/10))^2)/25 + (pi^2*abs(sin((timesteps(n)*pi)/10))^2)/25 + 1/4)^(1/2))];
%         vdest(:,n) = [-(pi*sin((pi*timesteps(n))/10))/5; (pi*cos((pi*timesteps(n))/10))/5; 0.5];
%         adest(:,n) = [-(pi^2*cos((pi*timesteps(n))/10))/50; -(pi^2*sin((pi*timesteps(n))/10))/50; 0];
%     end

end

function [disx, disW] = dis_gen(n)
    dis_gain_x = 20;
    dis_gain_w = 6;

    step_disturbance = 0;

    scenario1 = [0.13*atan(n/2); 0.2*atan(n/2); 0.26*atan(n/2)];
    scenario2 = [0.13; 0.2; 0.26];
    scenario3 = [0.13*cos(0.1*pi*n); 0.2*cos(0.1*pi*n); 0.26*cos(0.1*pi*n)];
    dis_fn = scenario2;
    
    disx = dis_gain_x * dis_fn;
    disW = dis_gain_w * dis_fn;

    if and(or(n < 5, n > 7), step_disturbance)
        disx = 0;
        disW = 0;
    end
end