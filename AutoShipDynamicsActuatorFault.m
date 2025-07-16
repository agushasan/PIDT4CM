%% Research code by Agus Hasan

clear;
clc;

%% Simulation Time
tf  = 20;
dt  = 0.001;
t   = dt:dt:tf;

%% System description
m    = 23.8;
Izz  = 1.76;
xg   = 0.046;
Xud  = -2;
Yvd  = -10;
Yrd  = 0;
Nvd  = 0;
Nrd  = -1;

%% Number of states
n = 6;

%% System parameters
M    = [m-Xud 0 0;0 m-Yvd m*xg-Yrd;0 m*xg-Nvd Izz-Nrd];

A   = eye(6);
B_t = [1 0 0;0 1 0; 0 0 1];
B   = dt*[0 0 0;0 0 0;0 0 0;inv(M)*B_t];
C   = eye(6);

%% noise
QF = 0.01*eye(rank(A));
RF = 1*eye(rank(C));

%% Initialization
x        = [4;0;0;0;0;0];
xhat     = [2;1;0;0;0;0];
xbar     = [1;-1;0;0;0;0];
Pplus    = eye(rank(A));
theta    = [0 0 0]';
thetahat = [0 0 0]';
thetabar = [0 0 0]';
temp     = [0 0 0 0 0 0]';
 
%% Paramater
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);

Xu   = -0.7225;
Xuu  = -1.3274;
Yv   = -0.8612;
Yvv  = -36.2823;
Yr   = 0.1079;
Nv   = 0.1052;
Nr  = -0.5;
Nrr = -1;

%% Control input
u     = [50 10 1]';

%% Parameters for AEKF
Psi         = -B*diag(u);
S           = 0.1;
UpsilonPlus = 0*B;
lambda      = 0.9;
a           = 0.999;

%% Parameters for AO
lambdav = 0.95;
lambdat = 0.9;
Rv      = 0.001*eye(n);
Rt      = 0.001*eye(n);
Pv      = 0.001*eye(n);
Pt      = 0.001*eye(3);
Gamma1  = zeros(6,3);

%% For plotting
uArray          = [];
xArray          = [];
xhatArray       = [];
xbarArray       = [];
thetaArray      = [];
thetahatArray   = [];
thetabarArray   = [];

%%
% Simulation
for i=1:(tf/dt)
    
    %u     = [50 10*cos(i*dt) 10]';
    
    Psi   = -B*diag(u);
    
    % permanent fault
    if i>4600
        theta(1) = 0.634;
    end
    
    % transient fault
    if i>6400
        theta(2) = 0.8+(0.8/(20000-6400))*(i-20000);
    end
    
    % intermittent fault
    if i>2000
        theta(3) = 0.424;
    end
    if i>2400
        theta(3) = 0;
    end
    if i>4800
        theta(3) = 0.321;
    end        
    if i>5100
        theta(3) = 0;
    end
    if i>8400
        theta(3) = 0.743;
    end
    if i>9800
        theta(3) = 0;
    end
    if i>14200
        theta(3) = 0.563;
    end
    if i>14500
        theta(3) = 0;
    end
    if i>16300
        theta(3) = 0.213;
    end
    if i>19500
        theta(3) = 0;
    end
    
    uArray         = [uArray u];
    xArray         = [xArray x];
    xhatArray      = [xhatArray xhat];
    thetaArray     = [thetaArray theta];
    thetahatArray  = [thetahatArray thetahat];
    xbarArray      = [xbarArray xbar];
    thetabarArray  = [thetabarArray thetabar]; 
    
    c13 = -m22*x(5)-((m23+m32)/2)*x(6);
    c23 = m11*x(4);
    Cv  = [0 0 c13; 0 0 c23;-c13 -c23 0];
    Dv  = -[Xu+Xuu*abs(x(4)) 0 0;0 Yv+Yvv*abs(x(5)) Yr;0 Nv Nr+Nrr*abs(x(6))];
    
    x = A*x+dt*[cos(x(3))*x(4)-sin(x(3))*x(5);sin(x(3))*x(4)+cos(x(3))*x(5);x(6);-inv(M)*(Cv+Dv)*[x(4);x(5);x(6)]]+B*u+Psi*theta+QF*dt*randn(6,1);
    y = C*x+RF*dt*randn(6,1);
    
    % Estimation using Adaptive Extended Kalman filter
    FX    = A+dt*[0 0 -sin(xhat(3))*xhat(4)-cos(xhat(3))*xhat(5) cos(xhat(3)) -sin(xhat(3)) 0; 0 0 cos(xhat(3))*xhat(4)-sin(xhat(3))*xhat(5) sin(xhat(3)) cos(xhat(3)) 0; 0 0 0 0 0 1;
                   -inv(M)*[0 0 0 Xu+2*Xuu*abs(xhat(4)) -m22*xhat(6) -m22*xhat(5)-(m23+m32)*xhat(6);0 0 0 Yr*xhat(6)+m11*xhat(6) Yv+2*Yvv*abs(xhat(5)) Yr+m11*xhat(4); 0 0 0 m22*xhat(5)+((m23+m32)/2)*xhat(6)-m11*xhat(5) m22*xhat(4)+Nv-m11*xhat(4) ((m23+m32)/2)*xhat(4)+Nr+2*Nrr*abs(xhat(6))]];
 
    Pmin  = FX*Pplus*FX'+QF;
    Sigma = C*Pmin*C'+RF;
    KF    = Pmin*C'*inv(Sigma);
    Pplus = (eye(rank(A))-KF*C)*Pmin;
     
    ytilde = y-C*xhat;
    QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
    RF    = a*RF + (1-a)*(ytilde*ytilde'+C*Pmin*C');
 
    Upsilon = (eye(rank(A))-KF*C)*FX*UpsilonPlus+(eye(rank(A))-KF*C)*Psi;
    Omega   = C*FX*UpsilonPlus+C*Psi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma   = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
     
    thetahat  = thetahat + Gamma*(y-C*xhat);
    xhat      = A*xhat+dt*[cos(xhat(3))*xhat(4)-sin(xhat(3))*xhat(5);sin(xhat(3))*xhat(4)+cos(xhat(3))*xhat(5);xhat(6);-inv(M)*(Cv+Dv)*[xhat(4);xhat(5);xhat(6)]]+B*u+Psi*thetahat+KF*(y-C*xhat)+Upsilon*Gamma*(y-C*xhat);

    % Estimation using Adaptive Observer
    Kv = Pv*inv(Pv+Rv);
    Kt = Pt*Gamma1'*inv(Gamma1*Pt*Gamma1'+Rt);
    Gamma1 = (eye(n)-Kv)*Gamma1;

    xbar = xbar+(Kv+Gamma1*Kt)*(y-xbar);
    thetabar = thetabar-Kt*(y-xbar);

    c13b = -m22*xbar(5)-((m23+m32)/2)*xbar(6);
    c23b = m11*xbar(4);
    Cvb  = [0 0 c13; 0 0 c23;-c13 -c23 0];
    Dvb  = -[Xu+Xuu*abs(xbar(4)) 0 0;0 Yv+Yvv*abs(xbar(5)) Yr;0 Nv Nr+Nrr*abs(xbar(6))];    
    
    xbar      = A*xbar+dt*[cos(xbar(3))*xbar(4)-sin(xbar(3))*xbar(5);sin(xbar(3))*xbar(4)+cos(xbar(3))*xbar(5);xbar(6);-inv(M)*(Cvb+Dvb)*[xbar(4);xbar(5);xbar(6)]]+B*u+Psi*thetabar;
    thetabar = thetabar;
    Pv = (1/lambdav)*eye(n)*(eye(n)-Kv)*Pv*eye(n);
    Pt = (1/lambdat)*(eye(3)-Kt*Gamma1)*Pt;
    Gamma1 = eye(n)*Gamma1-Psi;

end

figure(1)
plot(xArray(1,:),xArray(2,:), 'r', 'LineWidth', 6)
hold on;
plot(xhatArray(1,:),xhatArray(2,:), 'g:', 'LineWidth', 6)
hold on;
plot(xbarArray(1,:),xbarArray(2,:), 'b:', 'LineWidth', 6)
hold on;
plot(xArray(1,1),xArray(2,1), 'O', 'LineWidth', 20)
hold on;
plot(xArray(1,end),xArray(2,end), 'O', 'LineWidth', 20)
xlabel('x','FontSize',48)
ylabel('y','FontSize',48)
grid on;
grid minor;
set(gca,'FontSize',36)
legend('true trajectory','estimated trajectory AEKF','estimated trajectory AO','start','end','FontSize',48);

% figure(2)
% plot(t,uArray(1,:), 'r', 'LineWidth', 6)
% hold on;
% plot(t,uArray(2,:), 'g:', 'LineWidth', 6)
% hold on;
% plot(t,uArray(3,:), 'b:', 'LineWidth', 6)
% xlabel('x','FontSize',24)
% ylabel('y','FontSize',24)
% grid on;
% grid minor;
% set(gca,'FontSize',36)
% legend('\tau_u','\tau_v','\tau_r','FontSize',24);

figure(2)
subplot(3,1,1)
plot(t,thetaArray(1,:), 'r', 'LineWidth', 6)
hold on;
plot(t,thetahatArray(1,:), 'g:', 'LineWidth', 6)
hold on;
plot(t,thetabarArray(1,:), 'b:', 'LineWidth', 6)
ylabel('\theta_u','FontSize',48)
grid on;
grid minor;
set(gca,'FontSize',36)
ylim([-0.05 1])
subplot(3,1,2)
plot(t,thetaArray(2,:), 'r', 'LineWidth', 6)
hold on;
plot(t,thetahatArray(2,:), 'g:', 'LineWidth', 6)
hold on;
plot(t,thetabarArray(2,:), 'b:', 'LineWidth', 6)
ylabel('\theta_v','FontSize',48)
grid on;
grid minor;
set(gca,'FontSize',36)
ylim([-0.05 1])
subplot(3,1,3)
plot(t,thetaArray(3,:), 'r', 'LineWidth', 6)
hold on;
plot(t,thetahatArray(3,:), 'g:', 'LineWidth', 6)
hold on;
plot(t,thetabarArray(3,:), 'b:', 'LineWidth', 6)
legend('true \theta','estimated \theta AEKF','estimated \theta AO','FontSize',48);
ylabel('\theta_r','FontSize',48)
xlabel('t (s)','FontSize',48)
grid on;
grid minor;
set(gca,'FontSize',36)
ylim([-0.05 1])

figure(3);
subplot(3,1,1)
plot(t,thetaArray(1,:)-thetahatArray(1,:), ':g', 'LineWidth', 6)
hold on;
plot(t,thetaArray(1,:)-thetabarArray(1,:), ':b', 'LineWidth', 6)
grid on;
grid minor;
ylabel('Error \theta_u','FontSize',48)
set(gca,'FontSize',36)
legend('Error AEKF','Error AO','FontSize',48);
ylim([-1 1])
subplot(3,1,2)
plot(t,thetaArray(2,:)-thetahatArray(2,:), ':g', 'LineWidth', 6)
hold on;
plot(t,thetaArray(2,:)-thetabarArray(2,:), ':b', 'LineWidth', 6)
grid on;
grid minor;
ylabel('Error \theta_v','FontSize',48)
set(gca,'FontSize',36)
ylim([-1 1])
subplot(3,1,3)
plot(t,thetaArray(3,:)-thetahatArray(3,:), ':g', 'LineWidth', 6)
hold on;
plot(t,thetaArray(3,:)-thetabarArray(3,:), ':b', 'LineWidth', 6)
grid on;
grid minor;
ylabel('Error \theta_r','FontSize',48)
xlabel('t (s)','FontSize',48)
set(gca,'FontSize',36)
ylim([-1 1])