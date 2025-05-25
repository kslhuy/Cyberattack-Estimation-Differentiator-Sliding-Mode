clear all
close all
clc

%% Parameters
Ts = 1/1000;
t = 0:Ts:20;
pas = 1e-4;
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
%% Auxiliary output (unchanged)
% Relative degree is equal to 3 (number of differentiation of the output to obtain the input for the first time
r = 3;
gamma = 3;

c1 = C(1,:);
Ca_1 = [c1;c1*A;c1*A^(gamma-1)];
c2 = C(2,:);
Ca_2 = [c2;c2*A];
Ca = [Ca_1 ; Ca_2];

%% Initial conditions (unchanged)
x_p0 = [10 10 0];
x_f0 = [0 0 0];
xf(:,1) = x_f0;
xp(:,1) = x_p0;
x(:,1) = [x_p0(1)-x_f0(1)-li_1,x_p0(2)-x_f0(2),0];
t0 = 0;
t(1) = t0;

ya11_h(1) = 4;
ya12_h(1) = 2;
ya13_h(1) = 1;
ya14_h(1) = -2;
ya21_h(1) = 0;
ya22_h(1) = 1;
ya23_h(1) = 6;

x_hat(:,1) = x_f0';
ya_h(:,1) = [ya11_h(1) ya12_h(1) ya13_h(1) ya14_h(1) ya21_h(1) ya22_h(1) ya23_h(1)];
nu_h(1) = -3;

K11 = 0.8;
K12 = 0.9;
K13 = 0.99;
K21 = 1;
K22 = 2;

v_10 = 30;
v_20 = 30;

%% Adaptive parameters
lambda_min = 0.1;  % Minimum to ensure positivity
lambda_max = 100;  % Maximum to prevent excessive growth
z_max = 10;        % Bound for integral terms
dead_zone = 1e-3;  % Dead zone for robustness

%% Initialize adaptive parameters
lambda11 = zeros(1,Nt);
lambda12 = zeros(1,Nt);
lambda13 = zeros(1,Nt);
lambda14 = zeros(1,Nt);
z14 = zeros(1,Nt);

lambda21 = zeros(1,Nt);
lambda22 = zeros(1,Nt);
lambda23 = zeros(1,Nt);
z23 = zeros(1,Nt);

% Initial values
lambda11(1) = 10;
lambda12(1) = 10;
lambda13(1) = 10;
lambda14(1) = 10;
z14(1) = 0;

lambda21(1) = 10;
lambda22(1) = 10;
lambda23(1) = 10;
z23(1) = 0;

f_c_hat_new(1) = 0;

%% Main simulation loop
for i=1:Nt-1
    %% Simulation of the preceding vehicle (unchanged)
    ai_1(i) = 3*sin(2*pi/10*t(i));
    dxp(:,i) = Ai_1*xp(:,i)+Bi_1*ai_1(i);
    xp(:,i+1) = xp(:,i)+dxp(:,i)*(t(i+1)-t(i));

    %% Simulation of the following vehicle (unchanged)
    if t(i)<4
        nu(i) = 0;
    elseif t(i)>=4 && t(i)<8
        nu(i) = 4*sin(3*pi/2*t(i));
    elseif t(i)>=9 && t(i)<12
        nu(i) = -2;
    elseif t(i)>=12 && t(i)<16
        nu(i) = 5;
    else
        nu(i) = 0;
    end
    mu(i) = ai_1(i)+nu(i);
    usyn(i) = (h*k1*k2)*mu(i) - k2*xf(2,i) - (k1+h*k1*k2)*ya13_h(i)...
        +1/h*(1-h*k1*k2)*ya12_h(i) + k2/h*ya11_h(i)...
        - k2/h*(Li-li_1) - (k1+h*k1*k2)*nu_h(i);
    dxf(:,i) = Ai*xf(:,i)+Bi*sat(usyn(i),100);
    xf(:,i+1) = xf(:,i)+dxf(:,i)*(t(i+1)-t(i));

    %% Error dynamics (unchanged)
    dx(:,i) = A*x(:,i)+B*xf(2,i)+F*mu(i)+Delt-W*nu(i);
    x(:,i+1) = x(:,i)+dx(:,i)*(t(i+1)-t(i));
    y(:,i) = C*x(:,i);
    y1(i) = y(1,i);
    y2(i) = y(2,i);

    %% Sliding mode differentiator
    g1(i) = c1*A*A*(B*xf(2,i)+F*mu(i)+Delt);
    g2(i) = c1*A*(B*xf(2,i)+F*mu(i)+Delt);

    gamma1 = 3;
    gamma2 = 2;
    % d = 1e-5;
    d = 1e-4;
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

    %% Second output
    w20(i) = ya21_h(i)-y2(i);
    w21(i) = lambda21(i)*abs(w20(i))^((gamma2)/(gamma2+1))*sat(w20(i),d);
    dya21_h(i) = ya22_h(i)-w21(i)-K21*w20(i);
    ya21_h(i+1) = ya21_h(i)+dya21_h(i)*(t(i+1)-t(i));

    w22(i) = lambda22(i)*abs(w21(i))^((gamma2-1)/(gamma2))*sat(w21(i),d);
    dya22_h(i) = ya23_h(i)-w22(i)-K22*(ya22_h(i)-(ya22_h(i)-w21(i)-K21*w20(i)));
    ya22_h(i+1) = ya22_h(i)+dya22_h(i)*(t(i+1)-t(i));

    w23(i) = lambda23(i)*sat(w22(i),d);
    dya23_h(i) = -w23(i) + g2(i);
    ya23_h(i+1) = ya23_h(i)+dya23_h(i)*(t(i+1)-t(i));



    %% Sliding variables with dead zone
    s11 = deadzone(w10(i), dead_zone);
    s12 = deadzone(w11(i), dead_zone);
    s13 = deadzone(w12(i), dead_zone);
    s21 = deadzone(w20(i), dead_zone);
    s22 = deadzone(w21(i), dead_zone);

    %% Lambda update laws with bounds
    if i < Nt-1
        % First output
        lambda11_update = lambda11(i) + Ts * abs(s11)^(7/4) * sat(s11,d);
        lambda12_update = lambda12(i) + Ts * abs(s12)^(5/3) * sat(s12,d);
        lambda13_update = lambda13(i) + Ts * abs(s13)^(3/2) * sat(s13,d);
        z14_update = z14(i) + Ts * sign(s13);
        lambda14_update = lambda14(i) + Ts * z14(i) * s13;

        lambda11(i+1) = min(max(lambda11_update, lambda_min), lambda_max);
        lambda12(i+1) = min(max(lambda12_update, lambda_min), lambda_max);
        lambda13(i+1) = min(max(lambda13_update, lambda_min), lambda_max);
        z14(i+1) = min(max(z14_update, -z_max), z_max);
        lambda14(i+1) = min(max(lambda14_update, lambda_min), lambda_max);

        % Second output
        lambda21_update = lambda21(i) + Ts * abs(s21)^(5/3) * sat(s21,d);
        lambda22_update = lambda22(i) + Ts * abs(s22)^(3/2) * sat(s22,d);
        z23_update = z23(i) + Ts * sign(s22);
        lambda23_update = lambda23(i) + Ts * z23(i) * s22;

        lambda21(i+1) = min(max(lambda21_update, lambda_min), lambda_max);
        lambda22(i+1) = min(max(lambda22_update, lambda_min), lambda_max);
        z23(i+1) = min(max(z23_update, -z_max), z_max);
        lambda23(i+1) = min(max(lambda23_update, lambda_min), lambda_max);
    end

    %% Reconstruction of unknown input (unchanged)
    ksi_a_h1(i) = ya14_h(i)+g1(i);
    ksi_a_h2(i) = ya23_h(i)+g2(i);
    theta = [ksi_a_h1(i) ksi_a_h2(i)];
    x_hat(:,i+1) = inv(Ca'*Ca)*Ca'*[ya11_h(i);ya12_h(i);ya13_h(i);ya21_h(i);ya22_h(i)];

    H = (c1*(A^(gamma1-1))*W);
    nu_h(i+1) = pinv(H'*H)*H'*(-ksi_a_h1(i)+c1*A^gamma1*[ya11_h(i);ya12_h(i);ya13_h(i)]+g1(i));

    H2 = (c2*(A^(gamma2-1))*W);
    nu_h2(i+1) = pinv(H2'*H2)*H2'*(-ksi_a_h2(i)+c2*A^gamma2*[y1(i);ya21_h(i);ya22_h(i)]+g2(i));

    C_a_blabla = [(c1*(A^(gamma1))) ; (c2*(A^(gamma2)))];
    C_a_title = [(c1*(A^(gamma1-1))) ; (c2*(A^(gamma2-1)))];
    U = C_a_title * W;
    f_c_hat(i,:) = pinv(U'*U)*U'*(-theta + C_a_blabla*(x_hat(:,i)) + [g1(i) ; g2(i)] );

    theta_new = [ya12_h(i);ya13_h(i);ya14_h(i);ya22_h(i);ya23_h(i)];
    CaW = Ca*W;
    f_c_hat_new(i+1) = inv(CaW'*CaW)*CaW'*(-theta_new + Ca*A*(x_hat(:,i)) + Ca*(B*xf(2,i)+F*mu(i)+Delt));


    
end


plot_errors
plot_lamda
