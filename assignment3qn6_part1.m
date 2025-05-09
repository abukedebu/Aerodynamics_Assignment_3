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
total_mass_no_wing = 4.9007; %5.0628 %4.9007
N_rotor = 2;
g_mars = 3.73;
rho_mars = 0.02; %0.02 %0.0275 enough for no wing to achieve 90% power
T_per_prop = total_mass_no_wing*g_mars/N_rotor;
Omega_ingen = 2800 * 2 * pi/60;
R_ingen = 0.6;
gamma = 1.15;
Cdo = 0.02;
P_per_rotor = 330.87/N_rotor;
V_tip = Omega_ingen * R_ingen;
Omega = Omega_ingen * R_ingen /R;

A_body = 0.1;
C_D_body = 0.4;
wing_density = 74; %74
v_drone = 10;

C_T = 0.009399641320025; %T_per_prop/(1/2*rho_mars*pi*R^2*V_tip^2);
C_la = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');

n = 500;
n_2 = 100;
max_iter = 500;
tol = 1e-6;

wingspan_list = linspace(0, 3, n_2);
wing_chord_list = [0.2 0.4 0.6 0.8 1]; %[0.2 0.4 0.6 0.8 1]
beta = linspace(1e-5, pi, n);
alpha_wing = deg2rad(25); %initial angle of wing wrt body
D_body = 0.5 * rho_mars * C_D_body * v_drone^2 * A_body;

required_total_T = zeros(length(wingspan_list), length(wing_chord_list));
idrag_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
lift_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
fdrag_wing_list = zeros(length(wingspan_list), length(wing_chord_list));
phi_body = zeros(length(wingspan_list), length(wing_chord_list));
weight_list = zeros(length(wingspan_list), length(wing_chord_list));
for wing_chord_index = 1:length(wing_chord_list)
    wing_chord = wing_chord_list(wing_chord_index);

    for wingspan_index = 1:length(wingspan_list)
        wingspan = wingspan_list(wingspan_index);
    
        
        coeffs = zeros(n, n);
        
        for i_theta = 1:n
            theta = beta(i_theta);
            constant = -4 * wingspan / (C_la * wing_chord);
        
            for i_an = 1:n
                coeffs(i_theta, i_an) = constant * sin(i_an * theta) - (i_an * sin(i_an * theta) / sin(theta));
            end
        end
        
        % Solve for A_n coefficients
        b = (alpha_lo - alpha_wing) * ones(n, 1);
        A = coeffs \ b;
        A_n_alphas = A';
        
        for i_theta = 1:n
            alpha_induced(i_theta) = sum((1:n) .* A_n_alphas .* sin((1:n) .* beta)./sin(beta));
        end
        
        alpha_eff = alpha_wing - alpha_induced;
        %fprintf("min alpha: %f\n", rad2deg(min(alpha_eff)))
        %fprintf("max alpha: %f\n", rad2deg(max(alpha_eff)))
        C_Df = interp1(alpha_data, Cd_data, alpha_eff);
        y_bar = wingspan*(1-cos(beta))/2;
        
        C_l_wing = 4*wingspan/wing_chord * sum(A_n_alphas .* sin((1:n) .* beta));
        C_di_wing = C_l_wing .* alpha_induced;
        lift_wing = 0.5 * C_l_wing * rho_mars * wingspan * wing_chord * v_drone^2;
        idrag_wing = trapz(y_bar, 0.5 * C_di_wing * rho_mars * wing_chord * v_drone^2);
        fdrag_wing = trapz(y_bar, 0.5 * C_Df * rho_mars * wing_chord * v_drone^2);
        
        lift_wing_list(wingspan_index, wing_chord_index) = lift_wing;
        idrag_wing_list(wingspan_index, wing_chord_index) = idrag_wing;
        fdrag_wing_list(wingspan_index, wing_chord_index) = fdrag_wing;
        volume_wing = area_airfoil/1*wing_chord * wingspan;
        mass_wing = wing_density * volume_wing;
        total_mass = total_mass_no_wing + mass_wing;
        mass_list(wingspan_index, wing_chord_index) = total_mass;
        weight = total_mass * g_mars;
        weight_list(wingspan_index, wing_chord_index) = weight;

        vert_force = weight-lift_wing;
        vert_force_list(wingspan_index, wing_chord_index) = vert_force;
        hori_force = idrag_wing+fdrag_wing+D_body;
        hori_force_list(wingspan_index, wing_chord_index) = hori_force;
        phi_body(wingspan_index, wing_chord_index) = atan(hori_force/vert_force);
        required_total_T(wingspan_index, wing_chord_index) = sqrt((hori_force)^2 + (vert_force)^2);
    end
end

%% Power
required_T_per_rotor = required_total_T/N_rotor;
v_H = sqrt(required_T_per_rotor/(2*rho_mars*pi*R^2));
vi_result = ones(length(wingspan_list), length(wing_chord_list));
vi_result_old = vi_result;
for iter = 1:max_iter
    vi_result = v_H.^2 ./(sqrt((v_drone .* cos(phi_body)).^2 + (v_drone .* sin(phi_body)+vi_result).^2));
    if abs(vi_result_old-vi_result) < tol
        break
    end
end

% for wing_chord_index = 1:length(wing_chord_list)
%     for wingspan_index = 1:length(wingspan_list)
%         temp = [1 2*v_drone*phi_body(wingspan_index, wing_chord_index) v_drone^2 0 -v_H(wingspan_index, wing_chord_index)^4];
%         roots_result = roots(temp);
%         vi_result(wingspan_index, wing_chord_index) = roots_result(abs(imag(roots_result)) < 1e-4 & real(roots_result) > 0);
%     end
% end
P_ideal_per_rotor = required_T_per_rotor.*(v_drone.*sin(phi_body)+vi_result);
%% Define theta
a = 2; % 4*a + 3*b cannot be lesser than zero
b = -2;
theta_tip = 12/(4*a+3*b) * (C_T/(sigma*C_la) + 1/4*sqrt(C_T) + alpha_lo/3);
theta_blade = a*theta_tip + b*theta_tip*y;

%% BEM
V_c_list = v_drone .* sin(phi_body);
omega_list = linspace(155, 300, n_2);
omega_to_move = zeros(wingspan_index,wing_chord_index);

for wingspan_index = 1:length(wingspan_list)
    for wing_chord_index = 1:length(wing_chord_list)
        V_c = V_c_list(wingspan_index, wing_chord_index);
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
            alpha_blade_list(omega_index) = max((alpha_blade));
            C_l = C_la*(alpha_blade - alpha_lo);
            C_d = interp1(alpha_data, Cd_data, alpha_blade, 'linear', 'extrap');
            %fprintf("omega: %f max alpha: %f\n", omega_test, max(rad2deg(alpha_blade)))
            
            dTBE_dr = 0.5 * N_blades * rho_mars * c .* V_rel .^ 2 .* (C_l .* cos(phi) - C_d .* sin(phi)); 
            dTBE_dr(isnan(dTBE_dr)) = 0;
            TBE(omega_index) = trapz(r_list, dTBE_dr);
            
        end
        if TBE(end) ~= 0
            omega_to_move(wingspan_index,wing_chord_index) = interp1(TBE, omega_list, required_T_per_rotor(wingspan_index,wing_chord_index));
        end
        % fprintf("alpha: %f\n", rad2deg(interp1(omega_list, alpha_blade_list, omega_to_move(wingspan_index,wing_chord_index))))

    end
end
%% New power calculation
P_new_per_rotor = zeros(length(wingspan_list), length(wing_chord_list));
for wing_chord_index = 1:length(wing_chord_list)
    P0 = 1/8 * rho_mars * N_blades * c * Cdo .* omega_to_move(:,wing_chord_index).^3 .* R^4;
    P_new_per_rotor(:,wing_chord_index) = gamma * P_ideal_per_rotor(:,wing_chord_index) + P0;
end

%% The total power consumption and rotor angle vs flight speed from 0 to 12 m/s
v_var_drone_list = linspace(0, 12, n_2);
chosen_chord = 0.2;
chosen_wingspan = 0;

%% Figures
figure(1); clf;
hold on
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, required_total_T(:, wingchord_index), 'LineWidth', 2, ...
        'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
    hold on
end
grid on;
xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Required Thrust $T$ [N]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'best', 'Interpreter', 'latex');

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
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, rad2deg(phi_body(:,wingchord_index)), 'LineWidth', 2, ...
        'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
    hold on
end
grid on;
xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Angle of Body [deg]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'best', 'Interpreter', 'latex');


%%
figure(8); clf;
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, 2*P_new_per_rotor(:,wingchord_index), 'LineWidth', 2, ...
        'DisplayName', ['$c_{wing}$ = ' num2str(wing_chord_list(wingchord_index)) ' m'], 'Color', colors(wingchord_index))
    hold on
end
grid on;

set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$P_{total}$ [W]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'northwest', 'Interpreter', 'latex');
yticks(0:250:1500);
% plot(xlim, 0.9.*[138.1385 138.1385],'black')
%plot(xlim, 0.9 * [P_per_rotor P_per_rotor],'black')


%%

figure(9); clf;
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, weight_list(:,wingchord_index), 'LineWidth', 2, ...
        'DisplayName', ['$c_{wing}$ = ' num2str(wing_chord_list(wingchord_index))  ' m'], 'Color', colors(wingchord_index))
    hold on
end
grid on
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$W_{total}$ [N]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize',20);


%%

figure(10); clf;
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, mass_list(:,wingchord_index), 'LineWidth', 2, ...
        'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
    hold on
end
grid on
xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('mass $M$ [kg]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'best', 'Interpreter', 'latex');
plot(xlim, [total_mass_no_wing total_mass_no_wing],'black')

figure(11); clf;
for wingchord_index = 1:length(wing_chord_list)
    plot(wingspan_list, omega_to_move(:,wingchord_index), 'LineWidth', 2, ...
        'DisplayName', ['chord length = ' num2str(wing_chord_list(wingchord_index))], 'Color', colors(wingchord_index))
    hold on
end
grid on
xlabel('Wingspan Length $b$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('omega $\Omega$ [s$^{-1}$]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'best', 'Interpreter', 'latex');
plot(xlim, [Omega Omega],'black')
