clc; clear;
close all;

%% Set default figure properties to light mode
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultAxesGridColor', [0.15, 0.15, 0.15]);
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultFigurePosition',  [80, 50, 750, 450]);
fontsize_general = 18;
colors = ["#1f77b4", "#d1495b", "#3c1642"];

%% Constants and parameters
E_per_battery = 10/6;      % Wh per battery
E_init = 20; % Wh
m_per_battery = 0.047;     % kg per battery
m_computer = 1;            % kg

R_ingenuity = 0.6;
rho_mars = 0.02;
g_mars = 3.73;

P_ingenuity_per_rotor = 88.62314; 
m_fus_ingeniuty = 0.3;
m_no_fus_ingeniuty = 0.9 + 0.07 + 0.28 + 0.25;

Omega_ingenuity = 2800 / 60 * 2 * pi;
Cdo = 0.02;
gamma = 1.15;

added_batteries = linspace(1, 175, 175);

%% Chosen configuration
chord = 0.139317;
r = 1;
N_blades = 2;
N_rotor = 2;

m_total = zeros(size(added_batteries));
m_motors_arr = zeros(size(added_batteries));
m_propellers_arr = zeros(size(added_batteries));
m_fus_arr = zeros(size(added_batteries));
m_battery_arr = zeros(size(added_batteries));
P_final_total = zeros(size(added_batteries));
flight_time_hover = zeros(size(added_batteries));

%% Loop over battery counts
for i = 1:length(added_batteries)
    Omega = Omega_ingenuity * R_ingenuity / r;
    sigma = N_blades * chord / (pi * r);
    weight_per_prop = 0.07/4 * N_blades * (r / R_ingenuity);
    m_propellers = weight_per_prop * N_rotor;

    P_new_per_rotor = P_ingenuity_per_rotor;
    tolerance = 1e-4;
    error = 1;
    max_iter = 500;
    iter = 0;
    m_0 = 0;

    while error > tolerance && iter < max_iter
        iter = iter + 1;

        m_motors = 0.25 / 2 * P_new_per_rotor / P_ingenuity_per_rotor * N_rotor;
        m_battery = m_per_battery * added_batteries(i);
        m_no_fus = m_computer + m_motors + m_propellers + m_battery + 0.5;
        m_fus = m_no_fus * m_fus_ingeniuty / m_no_fus_ingeniuty;
        m_total_i = m_no_fus + m_fus;

        T_per_rotor = m_total_i * g_mars / N_rotor;
        P_0 = (1/8) * rho_mars * sigma * pi * Cdo * Omega^3 * r^5;
        P_per_rotor = ((T_per_rotor)^(3/2) / sqrt(2 * rho_mars * pi * r^2)) * gamma + P_0;

        error = abs(m_0 - m_motors);
        m_0 = m_motors;
        P_new_per_rotor = P_per_rotor;
    end

    m_total(i) = m_total_i;
    m_motors_arr(i) = m_motors;
    m_propellers_arr(i) = m_propellers;
    m_fus_arr(i) = m_fus;
    m_battery_arr(i) = m_battery;
    P_final_total(i) = N_rotor * P_new_per_rotor;

    E_total = (E_init + added_batteries(i) * E_per_battery) * 3600;
    flight_time_hover(i) = E_total / P_final_total(i) / 60;
end

%% Optional Plot: Flight time vs battery count

[max_time, idx_max] = max(flight_time_hover);
optimal_batteries = added_batteries(idx_max);
m_additional_batteries = m_battery_arr(idx_max);


% Find the last index where battery mass is just under 2 kg
idx_under_2kg = find(m_battery_arr < 2, 1, 'last');
batteries_under_2kg = added_batteries(idx_under_2kg);
flight_time_under_2kg = flight_time_hover(idx_under_2kg);


figure(1), clf;
ontsize_general = 20;


plot(added_batteries, flight_time_hover, 'LineWidth', 2, 'Color',colors(1), 'DisplayName','Total hover flight time');
hold on
plot(optimal_batteries, max_time, 'LineStyle','none', 'Marker','pentagram', 'MarkerFaceColor',colors(2), 'MarkerSize',16, ...
    'DisplayName',  sprintf('$N_{b,opt}$ = %.f', optimal_batteries), 'MarkerEdgeColor','none');


plot(added_batteries(42), flight_time_hover(42), 'LineStyle','none', 'Marker','pentagram', 'MarkerFaceColor',colors(3), 'MarkerSize',16, ...
    'DisplayName', sprintf('$N_{b,  max}$ = %.f', 42) , 'MarkerEdgeColor','none');
xline(batteries_under_2kg, 'Color',colors(3), 'LineStyle','--', 'LineWidth',2, 'DisplayName','Payload threshold')

%plot( )
xlim([0, 175]);
xticks(25:25:175);
ylim([0, 22]);
grid on;
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$N_{b}$', 'FontSize', 22, 'Interpreter','latex');
ylabel('$t$ [min]', 'FontSize', 22, 'Interpreter','latex');
legend('Location','southeast', Interpreter='latex')
%title('Hover Time vs. Battery Count (Two-Rotor, N_{blade} = 2, R = 1 m)', 'FontSize', 16);





fprintf('Optimum batteries: %d batteries \n', optimal_batteries);
fprintf('Optimum battery mass: %.3f kg\n', m_additional_batteries);
fprintf('Optimum hover flight time: %.2f minutes\n', max_time);


fprintf('\nMaximum batteries: %d batteries \n', batteries_under_2kg);
fprintf('Maximum battery mass: %.3f kg\n', m_battery_arr(idx_under_2kg));
fprintf('Maximum hover flight time: %.2f minutes\n', flight_time_under_2kg);