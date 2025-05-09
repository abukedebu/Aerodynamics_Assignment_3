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
fontsize_general = 22;

colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];

%% Load data
data = load('ishii_data_0.1_transition.txt');
alpha_lo = deg2rad(-1.223744);
Cl_data = data(:,2);
alpha_data = deg2rad(data(:,1));
Cd_data = data(:,3);
C_l_gradient = gradient(Cl_data, alpha_data);
C_la = interp1(alpha_data, C_l_gradient, alpha_lo);
area_airfoil = 0.04267605; % obtained from integration of area of airfoil

%% Parameters
N_blades = 2;
R = 1;
r_list = linspace(0.04, R - 0.01, 100);
y = r_list/R;

c = 0.139317; %0.139317
sigma = (N_blades .* c) ./ (pi * R);
M0 = gradient(Cl_data, alpha_data);
total_mass = 4.9007; %5.0628 %4.9007
N_rotor = 2;
g_mars = 3.73;
rho_mars = 0.02; %0.02 %0.0275 enough for no wing to achieve 90% power
T_per_prop = total_mass*g_mars/N_rotor;
Omega_ingen = 2800 * 2 * pi/60;
R_ingen = 0.6;
gamma = 1.15;
Cdo = 0.02;
P_per_rotor = 330.87/N_rotor;
V_tip = Omega_ingen * R_ingen;
Omega = Omega_ingen * R_ingen /R;

A_body = 0.1;
C_D_body = 0.4;
wing_density = 1e-1; %74

C_T = 0.009399641320025; %T_per_prop/(1/2*rho_mars*pi*R^2*V_tip^2);
C_la = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');

n = 500;
n_2 = 100;
max_iter = 500;
tol = 1e-6;

v_drone_list = linspace(0, 12, n_2);
wingspan_list = [0];
wing_chord_list = [1]; %[0.2 0.4 0.6 0.8 1]
beta = linspace(1e-5, pi, n);
alpha_wing = deg2rad(25); %initial angle of wing wrt body
D_body = 0.5 * rho_mars * C_D_body * v_drone_list.^2 * A_body;

required_total_T = zeros(length(wingspan_list), length(wing_chord_list));
idrag_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
lift_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
fdrag_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
phi_body = zeros(length(wingspan_list), length(wing_chord_list));
weight_list = zeros(length(wingspan_list), length(wing_chord_list));

weight = total_mass * g_mars;

vert_force = weight;
hori_force = D_body;
phi_body = atan(hori_force/vert_force);
required_total_T = sqrt((hori_force).^2 + (vert_force).^2);

%% Power
required_T_per_rotor = required_total_T/N_rotor;
v_H = sqrt(required_T_per_rotor/(2*rho_mars*pi*R^2));
vi_result = ones(size(v_H));
vi_result_old = vi_result;
% for iter = 1:max_iter
%     vi_result = v_H.^2 ./ sqrt((v_drone_list .* cos(phi_body)).^2 + (v_drone_list .* sin(phi_body) + vi_result).^2);
%     if abs(vi_result_old-vi_result) < tol
%         break
%     end
% end


for wingspan_index = 1:length(v_drone_list)
    temp = [1 2*v_drone_list(wingspan_index)*phi_body(wingspan_index) v_drone_list(wingspan_index)^2 0 -v_H(wingspan_index)^4];
    roots_result = roots(temp);
    vi_result(wingspan_index) = roots_result(abs(imag(roots_result)) < 1e-4 & real(roots_result) > 0);
end
P_ideal_per_rotor = required_T_per_rotor.*(v_drone_list.*sin(phi_body)+vi_result);
%% Define theta
a = 2; % 4*a + 3*b cannot be lesser than zero
b = -2;
theta_tip = 12/(4*a+3*b) * (C_T/(sigma*C_la) + 1/4*sqrt(C_T) + alpha_lo/3);
theta_blade = a*theta_tip + b*theta_tip*y;

%% BEM
V_c_list = v_drone_list .* sin(phi_body);
omega_list = linspace(155, 300, n_2);
omega_to_move = zeros(size(V_c_list));

for v_drone_index = 1:length(v_drone_list)
    V_c = V_c_list(v_drone_index);
    TBE = zeros(size(omega_list));
    lambda_i = 0.01 * ones(size(r_list)); %set to random value, during initialization
    lambda_c = V_c./(omega_list * R);
    F = ones(size(r_list));
    phi = ones(size(r_list));
    for iter = 1:max_iter
        phi_old = phi;
        phi = lambda_i ./ y;
        f = (N_blades / 2) * ((R - r_list) ./ (r_list .* sin(phi)));
        F = (2/pi) * acos(exp(-f));

        dCTmom = 8*(lambda_c + lambda_i) .* lambda_i .* y;
        lambda_i = -0.5 * (lambda_c + sigma .* C_la./(8*F) - sqrt((lambda_c.^2 + sigma .* C_la./(8*F)).^2 + sigma .* C_la .* (theta_blade - alpha_lo) .* y./(2*F)));
        %lambda_i = (sigma .* C_la) ./ (16 .* F) .* ...
                   % (sqrt(1 + (32 * F) ./ (sigma .* C_la) .* (theta_blade - alpha_lo) .* y) - 1);
        if max(abs(phi - phi_old)) < tol
            break;
        end
    end
    for omega_index = 1:length(omega_list)
        omega_test = omega_list(omega_index);
        phi = real(phi);
        V_rel = omega_test .* r_list ./ cos(phi);
        alpha_blade = real(theta_blade - phi);
        C_l = C_la*(alpha_blade - alpha_lo);
        C_d = interp1(alpha_data, Cd_data, alpha_blade, 'linear', 'extrap');
        % fprintf("omega: %f max alpha: %f\n", omega_test, max(rad2deg(alpha_blade)))

        dTBE_dr = 0.5 * N_blades * rho_mars * c .* V_rel .^ 2 .* (C_l .* cos(phi) - C_d .* sin(phi)); 
        dTBE_dr(isnan(dTBE_dr)) = 0;
        TBE(omega_index) = trapz(r_list, dTBE_dr);

    end
    if TBE(end) ~= 0
        omega_to_move(v_drone_index) = interp1(TBE, omega_list, required_T_per_rotor(v_drone_index));
    end
    % if v_drone_index < 5
    %     disp(TBE(10:20))
    %     disp(required_T_per_rotor(v_drone_index))
    %     disp(omega_to_move(v_drone_index))
    % end

    % fprintf("alpha: %f\n", rad2deg(interp1(omega_list, alpha_blade_list, omega_to_move(wingspan_index,wing_chord_index))))

end

%% New power calculation
P_new_per_rotor = zeros(length(wingspan_list), length(wing_chord_list));
P0 = 1/8 * rho_mars * N_blades * c * Cdo .* omega_to_move.^3 .* R^4;
P_new_per_rotor = gamma * P_ideal_per_rotor + P0;

%% Figures
wingchord_index = 1;

figure(1); clf;
hold on
plot(v_drone_list, required_total_T, 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
grid on;
xlabel('Drone Velocity $v$ [m/s]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Required Thrust $T$ [N]', 'Interpreter', 'latex', 'FontSize', 22);

% figure(2); clf;
% for wingchord_index = 1:length(wing_chord_list)
%     plot(wingspan_list, idrag_wing_list(:,wingchord_index), 'LineWidth', 2, ...
%         'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
%     hold on
% end
% grid on;
% xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('Induced Drag [-]', 'Interpreter', 'latex', 'FontSize', 22);
% legend('Location', 'best', 'Interpreter', 'latex');

% figure(3); clf;
% for wingchord_index = 1:length(wing_chord_list)
%     plot(wingspan_list, lift_wing_list(:,wingchord_index), 'LineWidth', 2, ...
%         'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
%     hold on
% end
% grid on;
% xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('Lift from Wing [N]', 'Interpreter', 'latex', 'FontSize', 22);
% legend('Location', 'best', 'Interpreter', 'latex');

% figure(4); clf;
% for wingchord_index = 1:length(wing_chord_list)
%     plot(wingspan_list, fdrag_wing_list(:,wingchord_index), 'LineWidth', 2, ...
%         'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
%     hold on
% end
% grid on;
% xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('Skin Friction Drag [N]', 'Interpreter', 'latex', 'FontSize', 22);
% legend('Location', 'best', 'Interpreter', 'latex');

% figure(5); clf;
% for wingchord_index = 1:length(wing_chord_list)
%     plot(wingspan_list, vert_force_list(:,wingchord_index), 'LineWidth', 2, ...
%         'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
%     hold on
% end
% grid on;
% xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('Vertical Force on Body [N]', 'Interpreter', 'latex', 'FontSize', 22);
% legend('Location', 'best', 'Interpreter', 'latex');

% figure(6); clf;
% for wingchord_index = 1:length(wing_chord_list)
%     plot(wingspan_list, hori_force_list(:,wingchord_index), 'LineWidth', 2, ...
%         'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
%     hold on
% end
% grid on;
% xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('Horizontal Force [N]', 'Interpreter', 'latex', 'FontSize', 22);
% legend('Location', 'best', 'Interpreter', 'latex');

figure(7); clf;
plot(v_drone_list, rad2deg(phi_body), 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
hold on
grid on;
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlim([0,12]);
xlabel('$V_{\infty}$ [m/s]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\beta$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
%%
figure(8); clf;
% plot(v_drone_list.*cos(phi_body)./v_H, N_rotor*P_new_per_rotor./(required_T_per_rotor .* v_H), 'LineWidth', 2, ...
plot(v_drone_list, N_rotor*P_new_per_rotor, 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
hold on
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');


grid on;
xlabel('$V_{\infty}$ [m/s]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$P$ [W]', 'Interpreter', 'latex', 'FontSize', 22);
% plot(xlim, 2*0.9.*[180.1286 180.1286],'black')
xlim([0,12])
plot(xlim, 2*[P_per_rotor P_per_rotor],'black', 'LineStyle','--', 'LineWidth',2)
%%
figure(9); clf;
plot(v_drone_list, P_ideal_per_rotor, 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
hold on
grid on;
xlabel('Drone Velocity $v$ [m/s]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('P ideal per rotor []', 'Interpreter', 'latex', 'FontSize', 22);
% plot(xlim, [P_per_rotor P_per_rotor],'black')

figure(10); clf;
plot(v_drone_list, omega_to_move, 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
hold on
grid on;
xlabel('Drone Velocity $v$ [m/s]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Omega [1/s]', 'Interpreter', 'latex', 'FontSize', 22);

figure(11); clf;
plot(v_drone_list.*cos(phi_body)./v_H, P_ideal_per_rotor./(required_T_per_rotor .* v_H), 'LineWidth', 2, ...
    'Color', colors(wingchord_index))
hold on
grid on;
xlabel('$V_{\inf}\cos(\beta)+v_i$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$P/Tv_H$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
% plot(xlim, [P_per_rotor P_per_rotor],'black')
