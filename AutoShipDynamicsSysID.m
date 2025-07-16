%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 20;
dt  = 0.001;
t   = dt:dt:tf;

n = 3;
r = 16;

%% True parameters

m    = 23.8;
Iz   = 1.76;
xg   = 0.046;

Xud  = -2;
Yvd  = -10;
Yrd  = 0;
Nvd  = 0;
Nrd  = -1;

Xu   = -0.7225;
Xuu  = -1.3274;
Yv   = -0.8612;
Yvv  = -36.2823;
Yr   = 0.1079;
Nv   = 0.1052;
Nr  = -0.5;
Nrr = -1;

%% system description
M = [m-Xud 0 0;0 m-Yvd m*xg-Yrd;0 m*xg-Nvd Iz-Nrd];
B = dt*[zeros(3);inv(M)];
C = [0 0 0 1 0 0;0 0 0 0 1 0; 0 0 0 0 0 1];

%% noise
R = 0.0001;

%% state initialization
x        = [3;3;0;0;0;0];
y        = [0;0;0];
vbar     = [0;0;0];
vhat     = [0;0;0];
thetabar = zeros(r,1);
thetahat = zeros(r,1);
 
%% known paramaters
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);
mt  = m22*m33-m23*m32;

alpha1 = dt*(m-Yvd)/(m-Xud);
alpha2 = dt*(m*xg-Yrd)/(m-Xud);
alpha3 = (-dt*(Iz-Nrd)*(m-Xud)/mt)+(dt*(m*xg-Yrd)*(m*xg-Yrd)/mt);
alpha4 = dt*(m*xg-Yrd)*(Xud-Yvd)/mt;
alpha5 = (dt*(m*xg-Nvd)*(m-Xud)/mt)-(dt*(m-Yvd)*(m*xg-Yrd)/mt);
alpha6 = -dt*(m-Yvd)*(Xud-Yvd)/mt;

beta1  = dt*Xu/(m-Xud);
beta2  = dt*Xuu/(m-Xud);
beta3  = (dt*(Iz-Nrd)*Yv/mt)-(dt*(m*xg-Yrd)*Nv/mt);
beta4  = (dt*(Iz-Nrd)*Yr/mt)-(dt*(m*xg-Yrd)*Nr/mt);
beta5  = dt*(Iz-Nrd)*Yvv/mt;
beta6  = -dt*(m*xg-Yrd)*Nrr/mt;
beta7  = (dt*(m-Yvd)*Nv/mt)-(dt*(m*xg-Nvd)*Yv/mt);
beta8  = (dt*(m-Yvd)*Nr/mt)-(dt*(m*xg-Nvd)*Yr/mt);
beta9  = -dt*(m*xg-Nvd)*Yvv/mt;
beta10 = dt*(m-Yvd)*Nrr/mt;

%% initial control inputs
u     = [40 10 1]';
%u     = [5 10 0]';

%% for plotting
uArray          = [];
xArray          = [];
yArray          = [];
vbarArray       = [];
vhatArray       = [];
thetabarArray   = [];
thetahatArray   = [];

%% Initialization for estimator

lambdav = 0.99;
lambdat = 0.9999;
Rv      = 0.001*eye(n);
Rt      = 0.001*eye(n);
Pv      = 0.001*eye(n);
Pt      = 0.001*eye(r);
Gamma   = zeros(n,r);

Pplus       = 10000000*eye(n);
QF          = 0.0001*eye(n);
RF          = 100000*eye(n);
a           = 0.999;
UpsilonPlus = 0*zeros(n,r);
S           = 1*eye(r);
lambda      = 0.999999;

%% simulation
for i=1:(tf/dt)
    
    u     = [40*cos(i*dt) 10*sin(i*dt) 1*sin(i*dt)]';

    uArray         = [uArray u];
    xArray         = [xArray x];
    yArray         = [yArray y];
    vbarArray      = [vbarArray vbar];
    vhatArray      = [vhatArray vhat];
    thetabarArray  = [thetabarArray thetabar]; 
    thetahatArray  = [thetahatArray thetahat]; 

    Cvv = [alpha1*x(5)*x(6)+alpha2*x(6)^2;alpha3*x(4)*x(6)+alpha4*x(4)*x(5);alpha5*x(4)*x(6)+alpha6*x(4)*x(5)];
    Dvv = [beta1*x(4)+beta2*abs(x(4))*x(4);beta3*x(5)+beta4*x(6)+beta5*abs(x(5))*x(5)+beta6*abs(x(6))*x(6);beta7*x(5)+beta8*x(6)+beta9*abs(x(5))*x(5)+beta10*abs(x(6))*x(6)];

    x = x+[dt*(cos(x(3))*x(4)-sin(x(3))*x(5));dt*(sin(x(3))*x(4)+cos(x(3))*x(5));dt*x(6);Cvv+Dvv]+B*u;
    y = C*x+R*rands(3,1);

    Phi = [y(2)*y(3) y(3)^2 y(1) abs(y(1))*y(1) 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 y(1)*y(3) y(1)*y(2) y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3) 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 y(1)*y(3) y(1)*y(2) y(2) y(3) abs(y(2))*y(2) abs(y(3))*y(3)];
    
    % Estimation using adaptive observer
    Kv = Pv*inv(Pv+Rv);
    Kt = Pt*Gamma'*inv(Gamma*Pt*Gamma'+Rt);
    Gamma = (eye(n)-Kv)*Gamma;

    vbar = vbar+(Kv+Gamma*Kt)*(y-vbar);
    thetabar = thetabar-Kt*(y-vbar);

    vbar = eye(n)*vbar+dt*inv(M)*(u)+Phi*thetabar;
    thetabar = thetabar;
    Pv = (1/lambdav)*eye(n)*(eye(n)-Kv)*Pv*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

    % Estimation using adaptive KF
    Pmin  = Pplus+QF;
    Sigma = Pmin+RF;
    KF    = Pmin*inv(Sigma);
    Pplus = (eye(n)-KF)*Pmin;
     
    ytilde = y-(vhat+dt*inv(M)*(u)+Phi*thetahat);
    QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
    RF    = a*RF + (1-a)*(ytilde*ytilde'+Pmin);
 
    Upsilon = (eye(n)-KF)*(eye(n))*UpsilonPlus+(eye(n)-KF)*Phi;
    Omega   = eye(n)*UpsilonPlus+Phi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma1  = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma1*ytilde;
    vhat      = vhat+dt*inv(M)*(u)+Phi*thetahat+KF*ytilde+Upsilon*Gamma1*ytilde;

end

Temp1bar = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(7,:);thetabarArray(13,:)];
Temp2bar = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetabarArray(8,:);thetabarArray(14,:)];
Temp1hat = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetahatArray(7,:);thetahatArray(13,:)];
Temp2hat = inv([-dt*m23/mt dt*m33/mt;dt*m22/mt -dt*m32/mt])*[thetahatArray(8,:);thetahatArray(14,:)];

% figure(1)
% plot(xArray(1,:),xArray(2,:), ':b', 'LineWidth', 6)
% hold on;
% plot(xArray(1,1),xArray(2,1), 'ok', 'LineWidth', 20)
% hold on;
% plot(xArray(1,end),xArray(2,end), 'or', 'LineWidth', 20)
% set(gca,'color','#EEEEEE','LineWidth',3,'FontSize',36)
% grid on;
% grid minor;
% legend('trajectory','start','end','FontSize',48)
% xlabel('x (m)','FontSize',48)
% ylabel('y (m)','FontSize',48)

figure(1)
subplot(3,1,1)
plot(t,yArray(1,:), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(1,1:100:end), 'g--', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vbarArray(1,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('u [m/s]','FontSize',48)
subplot(3,1,2)
plot(t,yArray(2,:), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(2,1:100:end), 'g--', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vbarArray(2,1:100:end), 'b:', 'LineWidth', 6)
grid on;
grid minor;
set(gca,'FontSize',36)
ylabel('v [m/s]','FontSize',48)
subplot(3,1,3)
plot(t,yArray(3,:), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vhatArray(3,1:100:end), 'g--', 'LineWidth', 6)
hold on;
plot(t(1:100:end),vbarArray(3,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('r [rad/s]','FontSize',48)
xlabel('time (s)','FontSize',48)
legend('measured','estimated AEKF','estimated AO','FontSize',48)

figure(2)
subplot(4,2,1)
plot(t,Xu*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetahatArray(3,1:100:end)/dt, 'g:', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetabarArray(3,1:100:end)/dt, 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylim([Xu-4 Xu+4]);
ylabel('X_u','FontSize',48)
subplot(4,2,2)
plot(t,Xuu*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetahatArray(4,1:100:end)/dt, ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),m11*thetabarArray(4,1:100:end)/dt, 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('X_{uu}','FontSize',48)
legend('true parameter','estimated AEKF','estimated AO','FontSize',36)
ylim([Xuu-4 Xuu+4]);
subplot(4,2,3)
plot(t,Nv*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1hat(1,1:100:end), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1bar(1,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylim([Nv-4 Nv+4]);
ylabel('N_v','FontSize',48)
subplot(4,2,4)
plot(t,Yv*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1hat(2,1:100:end), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp1bar(2,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('Y_v','FontSize',48)
ylim([Yv-4 Yv+4]);
subplot(4,2,5)
plot(t,Nr*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2hat(1,1:100:end), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2bar(1,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylim([Nr-4 Nr+4]);
ylabel('N_r','FontSize',48)
subplot(4,2,6)
plot(t,Yr*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2hat(2,1:100:end), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),Temp2bar(2,1:100:end), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('Y_r','FontSize',48)
ylim([Yr-4 Yr+4]);
subplot(4,2,7)
plot(t,Yvv*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),mt*thetahatArray(9,1:100:end)/(dt*m33), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),mt*thetabarArray(9,1:100:end)/(dt*m33), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
xlabel('time (s)','FontSize',48)
ylim([Yvv-4 Yvv+4]);
ylabel('Y_{vv}','FontSize',48)
subplot(4,2,8)
plot(t,Nrr*ones(length(t),1), 'r', 'LineWidth', 6)
hold on;
plot(t(1:100:end),-mt*thetahatArray(10,1:100:end)/(dt*m23), ':g', 'LineWidth', 6)
hold on;
plot(t(1:100:end),-mt*thetabarArray(10,1:100:end)/(dt*m23), 'b:', 'LineWidth', 6)
set(gca,'FontSize',36)
grid on;
grid minor;
ylabel('N_{rr}','FontSize',48)
xlabel('time (s)','FontSize',48)
ylim([Nrr-4 Nrr+4]);