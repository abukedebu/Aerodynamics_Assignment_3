clc;
clear;
close all;



%% Set default figure properties to light mode
close all;
set(0, 'DefaultFigureColor', 'w');     % Default figure background color (white)
set(0, 'DefaultAxesColor', 'w');   % Transparent axes background
set(0, 'DefaultAxesXColor', 'k');     % Default black X-axis
set(0, 'DefaultAxesYColor', 'k');     % Default black Y-axis
set(0, 'DefaultAxesGridColor', [0.15, 0.15, 0.15]); % Default grid line color
set(0, 'DefaultTextColor', 'k');      % Default text color (black)
set(0, 'DefaultFigurePosition',  [80, 50, 750, 450]);
fontsize_general = 20;
colors = ["#00a5cf", "#d1495b", "#3c1642", "#198754", "#1e5d8a", "#f6976d"];


%% Load data
data = load('ishii_data_0.1_transition.txt');
alpha_data = data(:,1) * pi / 180;
Cl_data = data(:,2);
Cd_data = data(:,3);
alpha_lo = -1.223744 * pi / 180;

%% Parameters
N_blades = 2;
R = 1;
r = linspace(0.04, R - 0.01, 100);
c = 0.139317; %0.139317
chord = c;
sigma = (N_blades .* chord) ./ (pi * R);
%alpha = 5 * pi/180;
alpha_L0 = alpha_lo;
M0 = gradient(Cl_data, alpha_data);
total_mass = 4.9007; %5.0628 %4.9007
N_rotor = 2;
g_mars = 3.73;
rho_mars = 0.02;
T_per_prop = total_mass*g_mars/N_rotor;
Omega_ingen = 2800 * 2 * pi/60;
R_ingen = 0.6;
gamma = 1.15;
Cdo = 0.02;
P_per_rotor = 330.87/N_rotor;
V_tip = Omega_ingen * R_ingen;
Omega = Omega_ingen * R_ingen /R;

C_T = T_per_prop/(1/2*rho_mars*pi*R^2*V_tip^2);
C_la = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');

y = r/R;

%% Define theta

%{
theta_tip = 12/5*(C_T./(sigma*C_la) + 1/4 * sqrt(C_T) + 1/3 * alpha_L0/3);
theta = theta_tip.*(2 - y);
%}


a = 2; % 4*a + 3*b cannot be lesser than zero
b = -2;
theta_tip = 12/(4*a+3*b) * (C_T/(sigma*C_la) + 1/4*sqrt(C_T) + alpha_L0/3);
theta = a*theta_tip + b*theta_tip*y;


%{
n = 1; % Be careful when n > 1.5
theta_tip = (3-n)*(C_T./(sigma*C_la) + 1/4*sqrt(C_T) + 1/3*alpha_L0);
theta = theta_tip./y.^n;
%}

%{
theta = 3*(C_T/(sigma*C_la) + 1/4*sqrt(C_T) + alpha_L0);
%}

theta_deg = rad2deg(theta);
%% Initialize
lambda_i = 0.1 * ones(size(r));
F = ones(size(r));
phi = zeros(size(r));
max_iter = 100;
tol = 1e-6;

%% Iteration
for iter = 1:max_iter
    phi_old = phi;

    phi = lambda_i ./ y;
    
    f = (N_blades / 2) * ((R - r) ./ (r .* sin(phi)));
    F = (2/pi) * acos(exp(-f));
    
    lambda_i = (sigma .* C_la) ./ (16 .* F) .* ...
               (sqrt(1 + (32 * F) ./ (sigma .* C_la) .* (theta - alpha_L0) .* y) - 1);
    
    if max(abs(phi - phi_old)) < tol
        break;
    end
end

alpha = theta - phi;
alpha_deg = rad2deg(alpha);
C_l = C_la*(alpha - alpha_L0);

V_rel = Omega .* r ./ cos(phi);
C_d = interp1(alpha_data, Cd_data, alpha);
dTBE_dr = 0.5 * N_blades * rho_mars * c .* V_rel .^ 2 .* (C_l .* cos(phi) - C_d .* sin(phi)); 
dCT = dTBE_dr/(0.5 * rho_mars * pi * R^2 .* (Omega * R).^2);
dP_dr = 0.5 * N_blades * Omega * rho_mars * c .* V_rel.^2 .* (C_l .* sin(phi) + C_d .* cos(phi)) .* r;
TBE = trapz(y*R, dTBE_dr);
P_new = trapz(y*R, dP_dr);

%% Plot dC_T vs r/R
figure(1), clf;
plot(r / R, real(dCT), 'Color', colors(1), 'LineWidth', 2)
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$r/R$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$dC_T/dr$ [1/m]', 'Interpreter', 'latex', 'FontSize', 22);
%title('Local Thrust Coefficient Distribution', 'Interpreter', 'latex');
grid on;

%% Plot dP vs r/R

figure(2), clf;
plot(r / R, real(dP_dr), 'Color', colors(1), 'LineWidth', 2);
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$r/R$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$dP/dr$ [W/m]', 'Interpreter', 'latex', 'FontSize', 22);
%title('Local Power Coefficient Distribution', 'Interpreter', 'latex');
grid on;


%% Plot alpha vs r/R
figure(4), clf;
plot(r / R, rad2deg(real(alpha)), 'Color', colors(1), 'LineWidth', 2);
hold on
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

% plot(xlim, rad2deg([alpha_L0 alpha_L0]), 'black')
% plot(xlim, [8 8], 'black') %stall angle
xlabel('$r/R$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
%title('Angle of Attack', 'Interpreter', 'latex');
grid on;

%% Plot twist
figure(5), clf;
plot(r / R, rad2deg(theta), 'Color', colors(1), 'LineWidth', 2);
hold on
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$r/R$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
%title('Twist Distribution', 'Interpreter', 'latex');
grid on;

%% Plot chord
figure(6), clf;
hold on
%xlim([0., 1])
%plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
ylim([0, 0.5]);
plot(r / R, chord * ones(size(r/R)), 'Color', colors(1), 'LineWidth', 2);
hold on
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$r/R$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('c [m]', 'Interpreter', 'latex', 'FontSize', 22);
%title('Chord Distribution', 'Interpreter', 'latex');
grid on;

% %% Plot lambda_i vs r/R
% figure;
% plot(r / R, real(lambda_i), 'b-', 'LineWidth', 2);
% xlabel('$r/R$', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('$\lambda_i$', 'Interpreter', 'latex', 'FontSize', 16);
% title('Induced Velocity Distribution', 'Interpreter', 'latex');
% grid on;

fprintf('Your T is: %.6f\n', TBE);
fprintf('Required T is: %.6f\n', T_per_prop);

fprintf('Your Power per rotor is: %.6f W\n', P_new);
fprintf('Your Power per rotor before is: %.6f W\n', P_per_rotor);