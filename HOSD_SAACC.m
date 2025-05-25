clear all
close all
clc

%% Parameters
Ts = 1/100;
t = 0:Ts:20;
% pas = 1e-4;
Nt = length(t);
tau = 0.4;
li_1 = 5;
Li = 12;
h = 0.5;
k1 = -0.8;
k2 = 2.5;
k3 = 2;

%% Continuous-time model of the preceeding vehicle
Ai_1 = [0 1 0;0 0 0;0 0 0];
Bi_1 = [0;1;0];
Ci_1 = zeros(1,3);
%% Continuous-time model of the following vehicle(usync as input)
Ai = [0 1 0;0 0 1;0 0 -1/tau];
Bi = [0;0;1/tau];
Ci = zeros(1,3);
%% Continuous-time model of SA-ACC system
A = [0 1 0;0 0 1;-k2/(tau*h) -k3/(tau*h) (k1-k3)/tau];
B = [0;0;k2/tau];
F = [0;0;k3/tau];
Delt = [0;0;k2*(Li-li_1)/(tau*h)];
W = [0;0;(k3-k1)/tau];
C = [1 0 0;0 1 0];
% C = [1 0 0];
%% Verifying rank condition
disp(' With actual output')
if rank(C*F)<rank(F)
    disp('Unknown input observer does not exist')
else
    disp('Unknown input observer exists')
end
%% Auxiliary output
% Relative degree is equal to 3 (number of differentiation of the output to obtain the input for the first time
r = 3;
gamma = 3;
% Augmented auxiliary output matrix
c1 = C(1,:);
Ca_1 = [c1;c1*A;c1*A^(gamma-1)];

c2 = C(2,:);
Ca_2 = [c2;c2*A];

Ca = [Ca_1 ; Ca_2];

%% Verifying rank condition
disp('With auxialiary output')
if rank(Ca_1*F)<rank(F)
    disp('Unknown input observer does not exist')
else
    disp('Unknown input observer exists')
end
%% System simulation
x_p0 = [10 10 0]; %Preceding car
x_f0 = [0 0 0];   % Following car
xf(:,1) = x_f0;
xp(:,1) = x_p0;
x(:,1) = [x_p0(1)-x_f0(1)-li_1,x_p0(2)-x_f0(2),0];
t0 = 0;
t(1) = t0;

%% Sliding mode Differentiator initial condition
% Auxiliary output initial guess
ya11_h(1) = 4;
ya12_h(1) = 2;
ya13_h(1) = 1;
ya14_h(1) = -2;
ya21_h(1) = 0;
ya22_h(1) = 1;
ya23_h(1) = 6;

x_hat(:,1) = x_f0';

ya_h(:,1) = [ya11_h(1) ya12_h(1) ya13_h(1) ya14_h(1) ya21_h(1) ya22_h(1) ya23_h(1)];

% Unknown input initial guess
nu_h(1) = -3;

K11 = 0.8;
K12 = 0.9;
K13 = 0.99;


K21 = 1;
K22 = 2;


v_10 = 30;
v_20 = 30;
t01 = 0;

lambda11(1) = 10;
lambda12(1) = 10;
lambda13(1) = 10;
lambda14(1) = 90;


lambda21(1) = 10;
lambda22(1) = 30;
lambda23(1) = 100;

f_c_hat_new(1) = 0;

for i=1:Nt-1
    %% Simulation of the preceeding vehicle
    ai_1(i) = 3*sin(2*pi/10*t(i));
    dxp(:,i) = Ai_1*xp(:,i)+Bi_1*ai_1(i);
    xp(:,i+1) = xp(:,i)+dxp(:,i)*(t(i+1)-t(i));
    %% Simulation of the following vehicle
    % Fault : Attack
    if t(i)<4
        nu(i) = 0;
    elseif t(i)>=4 && t(i)<8
%         nu(i) = 3*sin(3*pi/2*t(i))*exp(-0.01*(t(i)));
        
        nu(i) = 4*sin(3*pi/2*t(i));

    elseif t(i)>=9 && t(i)<12
%         nu(i) =  exp(-2*sin(3*pi/2*t(i))*0.01*(t(i)));
        nu(i) =  -2;
   elseif t(i)>=12 && t(i)<16
%         nu(i) =  exp(-2*sin(3*pi/2*t(i))*0.01*(t(i)));
        nu(i) =  5;

    else
        nu(i) =  0;%3;
    end
    mu(i) = ai_1(i)+nu(i);
    % Control law
    usyn(i) = (h*k1*k2)*mu(i) - k2*xf(2,i) - (k1+h*k1*k2)*ya13_h(i)...
        +1/h*(1-h*k1*k2)*ya12_h(i) + k2/h*ya11_h(i)...
        - k2/h*(Li-li_1) - (k1+h*k1*k2)*nu_h(i); 
    
%     usyn(i) = (h*k1*k2)*mu(i) - k2*xf(2,i) - (k1+h*k1*k2)*x_hat(3,i)...
%         +1/h*(1-h*k1*k2)*x_hat(2,i) + k2/h*x_hat(1,i)...
%         - k2/h*(Li-li_1) - (k1+h*k1*k2)*f_c_hat_new(i); 
    
%     usyn(i) = (h*k1*k2)*mu(i) - k2*xf(2,i) - (k1+h*k1*k2)*x(3,i)...
%         +1/h*(1-h*k1*k2)*x(2,i) + k2/h*x(1,i)...
%         - k2/h*(Li-li_1) - (k1+h*k1*k2)*nu(i);  
    
    dxf(:,i) = Ai*xf(:,i)+Bi*sat(usyn(i),100);
    xf(:,i+1) = xf(:,i)+dxf(:,i)*(t(i+1)-t(i));
    %% Simulation of the closed loops dynamics : error dynamics
    % Dynamic equation
    dx(:,i) = A*x(:,i)+B*xf(2,i)+F*mu(i)+Delt-W*nu(i);
    x(:,i+1) = x(:,i)+dx(:,i)*(t(i+1)-t(i));
    % Output equation
    y(:,i) = C*x(:,i) ;
    y1(i) = y(1,i) ;
    y2(i) = y(2,i) ;
    
    
%     y1(i) = y(1,i) +  rand*0.1;
%     y2(i) = y(2,i) + rand*0.01;
%     
    %% Sliding mode differentiator
    
    g1(i) = c1*A*A*(B*xf(2,i)+F*mu(i)+Delt);
    g2(i) = c1*A*(B*xf(2,i)+F*mu(i)+Delt);

    
    gamma1 = 3; % 1st Relative degree
    gamma2 = 2; % 2nd Relative degree
    
    
    d = 1e-5;
    
    %% Differentiator gains:
    %------------------
    lambda11(i) = 10;
    lambda12(i) = 11;
    lambda13(i) = 12;
    lambda14(i) = 85;
    %------------------
    lambda21(i) = 10;
    lambda22(i) = 15;
    lambda23(i) = 90;

    %% First output
    w10(i) = ya11_h(i)-y1(i);
    w11(i) = lambda11(i)*abs(w10(i))^((gamma1)/(gamma1+1))*sat(w10(i),d);
    dya11_h(i) = ya12_h(i)-w11(i)-K11*w10(i);
    ya11_h(i+1) = ya11_h(i)+dya11_h(i)*(t(i+1)-t(i));
    
    w12(i) = lambda12(i)*abs(w11(i))^((gamma1-1)/(gamma1))*sat(w11(i),d);
    dya12_h(i) = ya13_h(i)-w12(i)-K12*(ya12_h(i)-(ya12_h(i)-w11(i)-K11*w10(i)));
    ya12_h(i+1) = ya12_h(i)+dya12_h(i)*(t(i+1)-t(i));
    
    w13(i) = lambda13(i)*abs(w12(i))^((gamma1-2)/(gamma1-1))*sat(w12(i),d);
    dya13_h(i) = ya14_h(i)-w13(i)+g1(i)-K13*(ya13_h(i)-(ya13_h(i)-w12(i)-K12*(ya12_h(i)-(ya12_h(i)-w11(i)-K11*w10(i)))));
    ya13_h(i+1)= ya13_h(i)+dya13_h(i)*(t(i+1)-t(i));
    
    w14(i) = lambda14(i)*sat(w13(i),d);
    dya14_h(i) = -w14(i);
    ya14_h(i+1) = ya14_h(i)+dya14_h(i)*(t(i+1)-t(i));
    
    %% Second Output (N'est pas nécessaire)
    
    w20(i) = ya21_h(i)-y2(i);
    w21(i) = lambda21(i)*abs(w20(i))^((gamma2)/(gamma2+1))*sat(w20(i),d);
    dya21_h(i) = ya22_h(i)-w21(i)-K21*w20(i);
    ya21_h(i+1) = ya21_h(i)+dya21_h(i)*(t(i+1)-t(i));
    
    w22(i) = lambda22(i)*abs(w21(i))^((gamma2-1)/(gamma2))*sat(w21(i),d);
    dya22_h(i) = ya23_h(i)-w22(i)- K22*(ya22_h(i)-(ya22_h(i)-w21(i)-K21*w20(i)));
    ya22_h(i+1) = ya22_h(i)+dya22_h(i)*(t(i+1)-t(i));
    
    w23(i) = lambda23(i)*sat(w22(i),d);
    dya23_h(i) = -w23(i) +  g2(i) ;
    ya23_h(i+1) = ya23_h(i)+dya23_h(i)*(t(i+1)-t(i));
    
    
    %% Reconstruction of unknown input
    ksi_a_h1(i) = ya14_h(i)+g1(i);
    ksi_a_h2(i) = ya23_h(i)+g2(i); % HUY
    
    theta = [ksi_a_h1(i) ksi_a_h2(i)];  % HUY

    x_hat(:,i+1) = inv(Ca'*Ca)*Ca'*[ya11_h(i);
            ya12_h(i);
            ya13_h(i);
            ya21_h(i);
            ya22_h(i)];
    
        
        
        
    H = (c1*(A^(gamma1-1))*W);
    nu_h(i+1) = pinv(H'*H)*H'*(-ksi_a_h1(i)+c1*A^gamma1*[ya11_h(i);ya12_h(i);ya13_h(i)]+g1(i));
  
    % HUY 2ere output

    H2 = (c2*(A^(gamma2-1))*W);

    nu_h2(i+1) = pinv(H2'*H2)*H2'*(-ksi_a_h2(i)+c2*A^gamma2*[y1(i);ya21_h(i);ya22_h(i)]+g2(i));

    % HUY Formul ZHU
    C_a_blabla = [(c1*(A^(gamma1))) ; (c2*(A^(gamma2)))];
    C_a_title = [(c1*(A^(gamma1-1))) ; (c2*(A^(gamma2-1)))];
    U = C_a_title * W;
%     U = Ca * W;
    
    f_c_hat(i,:) = pinv(U'*U)*U'*(-theta + C_a_blabla*(x_hat(:,i)) + [g1(i) ; g2(i)] );
    
    %%%%     % EXACT formule

    theta_new = [ya12_h(i);
                ya13_h(i);
                ya14_h(i);
                ya22_h(i);
                ya23_h(i)];
    
    CaW = Ca*W;
    f_c_hat_new(i+1) = inv(CaW'*CaW)*CaW'*(-theta_new + Ca*A*(x_hat(:,i)) + Ca*(B*xf(2,i)+F*mu(i)+Delt) );
    

    
    
    
end
%% Results
figure(1)
set(gcf,'color','w');
subplot(3,1,1)
plot(t,ya11_h,'r-.',t,x(1,:),'b--','LineWidth',1.5);
% hold on
% plot(t,x_hat(1,:),'r-.',t,x(1,:),'b--','LineWidth',1.5);


grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('$$\delta_i$$','Interpreter','LaTeX','Fontsize',12);
d1 = legend('$$\hat{\delta_i}$$','$$\delta_i$$');
d1.Interpreter = "latex";
d1.FontSize = 18;
d11 = title('Relative position between host and preceding vehicle and its estimation');
d11.Interpreter = "latex";
d11.FontSize = 10;


subplot(3,1,2)
set(gcf,'color','w');
% plot(t,ya12_h,'r-.',t,x(2,:),'b--','LineWidth',1.5);
plot(t,x_hat(2,:),'r',t,x(2,:),'b--','LineWidth',1.5);

grid on

xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('\boldmath$\dot{\delta}_i$','Interpreter','LaTeX','Fontsize',12);
d2 = legend('$\dot{\hat{\delta}_i}$','$$\boldmath{\dot{\delta}_i}$$');
d2.Interpreter = "latex";
d2.FontSize = 18;
d12 = title('Relative velocity between host and preceding vehicle and its estimation');
d12.Interpreter = "latex";
d12.FontSize = 10;


subplot(3,1,3)
set(gcf,'color','w');
% plot(t,ya13_h,'r-.',t,x(3,:),'b--','LineWidth',1.5);
plot(t,x_hat(3,:),'r',t,x(3,:),'b--','LineWidth',1.5);

grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('\boldmath$\ddot{\delta}_i$','Interpreter','LaTeX','Fontsize',12);
d3 = legend('$$\ddot{\hat{\delta}_i}$$','$$\boldmath{\ddot{\delta}_i}$$');
d3.Interpreter = "latex";
d3.FontSize = 18;
d13 = title('Relative acceleration between host and preceding vehicle and its estimation');
d13.Interpreter = "latex";
d13.FontSize = 10;


figure(2)
set(gcf,'color','w');
rgb1 = [0.4660 0.6740 0.1880];
rgb2 = [0 0.4470 0.7410];
plot(t(2:end),nu,'Color',rgb1,'LineWidth',3);grid on
hold on;
plot(t, nu_h ,'Color',rgb2,'LineWidth',3);

% hold on ; 
% plot(t(2:end), f_c_hat(:,1) ,'r','LineWidth',1.5);
% plot(t(2:end), f_c_hat(:,2) ,'b','LineWidth',1.5);
% 
% plot(t, f_c_hat_new ,'g','LineWidth',1.5);

% plot(t(2:end), f_c_hat_new - f_c_hat(:,1)  ,'y','LineWidth',1.5);


grid on
xlabel('{Time} ${[sec]}$','Interpreter','LaTeX','Fontsize',12);
ylabel('$$f_{c}$$','Interpreter','LaTeX','Fontsize',12);
ylim([-90 30 ]);
d4 = legend('${\hat{f}_c}$','$$f_{c}$$');
d4.Interpreter = "latex";
d4.FontSize = 18;
d14 = title('Cyber-attack and its estimation');
d14.Interpreter = "latex";
d14.FontSize = 10;



ax = axes;

set(ax,'units','normalized','position',[0.2,0.2,0.2,0.2])
box(ax,'on')
plot(t,nu_h,'Color',rgb1,'LineWidth',1.5);grid on
hold on;
plot(t,[nu nu_h(end)] ,'Color',rgb2,'LineWidth',1.5);grid on
set(ax,'xlim',[0,2],'ylim',[-85,20])


