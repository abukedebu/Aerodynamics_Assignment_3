clc; clear;
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
fontsize_general = 18;
colors = ["#00a5cf", "#d1495b", "#3c1642", "#198754", "#1e5d8a", "#f6976d"];
%% Constants and parameters
N_rotor = 2; % change to 4 for quadcopter
m_battery = 0.5; % kg
E_per_battery = 20; % Wh
m_payload = 2; % kg
m_computer = 1; % kg

R = linspace(0.1, 1.5, 100);
R_ingenuity = 0.6;

rho_mars = 0.02;
g_mars = 3.73;

P_ingenuity_per_rotor = 88.62314; %134 %250 % 88.62314; % W
m_fus_ingeniuty = 0.3;
m_no_fus_ingeniuty = 0.9 + 0.07 + 0.28 + 0.25;

Omega_ingenuity = 2800 / 60 * 2 * pi;
Cdo = 0.02;
gamma = 1.15;

blade_counts = [2, 3, 4];
rotor_counts = [2, 4];


%% Try different numbers of blades and quadcopter
P_final_total = zeros(length(R), length(blade_counts), length(rotor_counts));
m_final = zeros(length(R), length(blade_counts), length(rotor_counts));
m_motors_arr = zeros(length(R), length(blade_counts), length(rotor_counts));
m_propellers_arr = zeros(length(R), length(blade_counts), length(rotor_counts));
m_fus_arr = zeros(length(R), length(blade_counts), length(rotor_counts));

for pos_rotor = 1:length(rotor_counts)
    N_rotor = rotor_counts(pos_rotor);
    for pos_blade = 1:length(blade_counts)
        N_blades = blade_counts(pos_blade);
        chord = 0.139317;
        for i = 1:length(R)
            r = R(i);
            Omega = Omega_ingenuity * R_ingenuity / r;
            sigma = N_blades*chord/(pi*r);
            weight_per_prop(i) = 0.07/4 * N_blades * (r / R_ingenuity);
            m_propellers = weight_per_prop(i) * N_rotor;
    
            P_new_per_rotor = P_ingenuity_per_rotor;
            tolerance = 1e-4;
            error = 1;
            max_iter = 500;
            iter = 0;
            m_0 = 0;
    
            while error > tolerance && iter < max_iter
                iter = iter + 1;
    
                m_motors = 0.25 / 2 * P_new_per_rotor / P_ingenuity_per_rotor * N_rotor;
                m_no_fus = m_computer + m_motors + m_propellers + m_battery + m_payload;
                m_fus = m_no_fus * m_fus_ingeniuty / m_no_fus_ingeniuty;
                m_motors_arr(i, pos_blade, pos_rotor) = m_motors;
                m_propellers_arr(i, pos_blade, pos_rotor) = m_propellers;
                m_fus_arr(i, pos_blade, pos_rotor) = m_fus;
                m_total = m_no_fus + m_fus;
    
                T_per_rotor = m_total * g_mars/N_rotor;
                P_0 = 1/8 * rho_mars * sigma * pi * Cdo * Omega^3 * r^5;
                P_per_rotor = ( (T_per_rotor)^(3/2) / sqrt(2 * rho_mars * pi * r^2) * gamma + P_0);
    
                error = abs(m_0 - m_motors);
                P_new_per_rotor = P_per_rotor;
                m_0 = m_motors;
            end
            
            m_final(i, pos_blade, pos_rotor) = m_total;
            P_final_total(i, pos_blade, pos_rotor) = N_rotor * P_new_per_rotor;
        end
    
    end
end

% N_blade = 2, N_rotor = 2, R = 1 m
[~, idx_R1] = min(abs(R - 1));
idx_blade = 1;
idx_rotor = 1;

m_comp_R1 = m_computer;
m_motors_R1 = m_motors_arr(idx_R1, idx_blade, idx_rotor);
m_propellers_R1 = m_propellers_arr(idx_R1, idx_blade, idx_rotor);
m_batt_R1 = m_battery;
m_payl_R1 = m_payload;
m_fus_R1 = m_fus_arr(idx_R1, idx_blade, idx_rotor);
m_total_R1 = m_final(idx_R1, idx_blade, idx_rotor);
P_final_total_R1 = P_final_total(idx_R1, idx_blade, idx_rotor);

% N_blade = 4, N_rotor = 2, R = 0.65 m
[~, idx_R2] = min(abs(R - 0.65));
idx_blade = 1;
idx_rotor = 2;
m_comp_R2 = m_computer;
m_motors_R2 = m_motors_arr(idx_R2, idx_blade, idx_rotor);
m_propellers_R2 = m_propellers_arr(idx_R2, idx_blade, idx_rotor);
m_batt_R2 = m_battery;
m_payl_R2 = m_payload;
m_fus_R2 = m_fus_arr(idx_R2, idx_blade, idx_rotor);
m_total_R2 = m_final(idx_R2, idx_blade, idx_rotor);
P_final_total_R2 = P_final_total(idx_R2, idx_blade, idx_rotor);

masses_R1 = [m_comp_R1, m_motors_R1, m_propellers_R1, m_batt_R1, m_payl_R1, m_fus_R1];
masses_R2 = [m_comp_R2, m_motors_R2, m_propellers_R2, m_batt_R2, m_payl_R2, m_fus_R2];
labels_text = {'Computer', 'Motors', 'Propellers', 'Battery', 'Payload', 'Fuselage'};

%%
 
figure (3), clf;
labels = cell(1, length(masses_R1));
for i = 1:length(masses_R1)

    labels{i} = sprintf('%s: %.2f kg', labels_text{i}, masses_R1(i));
end
clc
pie1 = pie(masses_R1, labels);

for i = 1:6
    j = 2*i;
    %get(pie1(j), 'Type')
    pie1(j-1).FaceColor = colors(i);
    pie1(j).FontSize = 22;
    pie1(j).Color = colors(i); 
    pie1(j).Interpreter = 'latex';
end
title('Two-rotor design, $R$ = 1m, $N_{blades}$ = 2', 'FontSize', 22, 'Interpreter','latex');

%%
figure (4), clf;
labels = cell(1, length(masses_R2));
for i = 1:length(masses_R2)
    labels{i} = sprintf('%s: %.2f kg', labels_text{i}, masses_R2(i));
end
pie2 = pie(masses_R2, labels);
for i = 1:6
    j = 2*i;
    %get(pie1(j), 'Type')
    pie2(j-1).FaceColor = colors(i);
    pie2(j).FontSize = 22;
    pie2(j).Color = colors(i); 
    pie2(j).Interpreter = 'latex';

end
title('Quadcopter, $R$ = 0.65m, $N_{blades}$ = 2', 'FontSize', 22, 'Interpreter','latex');

%title('Mass Breakdown at R = 0.65 m (2 Blades, 4 Rotors)', 'FontSize', 16);

%% plotting figure of power vs r
figure(1), clf;
set(gcf, "Position", [80, 50, 950, 450])
hold on;
%helicopter, 2 propeller
for pos_blade = 1:length(blade_counts)
    plot(R, P_final_total(:, pos_blade, 1), 'LineWidth', 2, ...
        'DisplayName', ['$N_{\mathrm{blade}}$=', num2str(blade_counts(pos_blade)), ', 2 rotors'], ...
        'Color', colors(pos_blade));
end

%quadcopter
for pos_blade = 1:length(blade_counts)
    plot(R, P_final_total(:, pos_blade, 2), 'LineWidth', 2, 'LineStyle', '--', ...
        'DisplayName', ['$N_{\mathrm{blade}}$=', num2str(blade_counts(pos_blade)), ', quadcopter'], ...
        'Color', colors(pos_blade));
end

grid on;
xlabel('$R$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$P$ [W]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'eastoutside', 'Interpreter','latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([0.1, 1.5]);
ylim([200, 1300]);
%xticks(-1:0.2:1);
plot(xlim, [min(P_final_total(P_final_total == min(P_final_total))), min(P_final_total(P_final_total == min(P_final_total)))], 'LineWidth', .1, 'HandleVisibility', 'off', 'Color', 'k');

%% plotting figure mass vs payload
figure(2), clf;
set(gcf, "Position", [80, 50, 950, 450])
hold on;
%helicopter, 2 propeller
for pos_blade = 1:length(blade_counts)
    plot(R, m_final(:, pos_blade, 1), 'LineWidth', 2, ...
        'DisplayName', ['$N_{\mathrm{blade}}$=', num2str(blade_counts(pos_blade)), ', 2 rotors'], ...
        'Color', colors(pos_blade));
end

%quadcopter
for pos_blade = 1:length(blade_counts)
    plot(R, m_final(:, pos_blade, 2), 'LineWidth', 2, 'LineStyle', '--', ...
        'DisplayName', ['$N_{\mathrm{blade}}$=', num2str(blade_counts(pos_blade)), ', quadcopter'], ...
        'Color', colors(pos_blade));
end

grid on;
xlabel('$R$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$M_{tot}$ [kg]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'eastoutside', 'Interpreter','latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([0.1, 1.5]);
ylim([4.5, 8]);
%xticks(-1:0.2:1);
plot(xlim, [min(m_final(m_final == min(m_final))), min(m_final(m_final == min(m_final)))], 'LineWidth', .1, 'HandleVisibility', 'off', 'Color', 'k');
