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
fontsize_general = 20;
colors = ["#00a5cf", "#d1495b", "#3c1642", "#198754", "#1e5d8a", "#f6976d"];

%% ------------------ GEOMETRY --------------------------------

ishii_geometry = load("ishii_geometry.dat");
ishii_x = ishii_geometry(:, 1);
ishii_y = ishii_geometry(:, 2);

break_idx = round(length(ishii_geometry(:, 1)) /2);

ishii_upper = ishii_geometry(1:break_idx, :);
ishii_lower =  flip(ishii_geometry(break_idx: end, :), 1);
ishii_camber = [ishii_upper(:, 1) mean([ishii_lower(:, 2) ishii_upper(:, 2)], 2)];


ubd_geometry = load("UBD5494_Geometry.dat");
ubd_x = ubd_geometry(:, 1);
ubd_y = ubd_geometry(:, 2);

break_idx = round(length(ubd_geometry(:, 1)) /2);

ubd_upper = ubd_geometry(1:break_idx, :);
ubd_lower =  flip(ubd_geometry(break_idx: end, :), 1);
ubd_camber = [ubd_upper(:, 1) mean([ubd_lower(:, 2) ubd_upper(:, 2)], 2)];


%%
clc;
figure(1), clf;
hold on
xlim([0, 1])
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');

plot(ubd_x, ubd_y, 'Color', colors(2), 'LineStyle','-', 'LineWidth',2, 'DisplayName', 'UBD-5459');
plot(ishii_x, ishii_y, 'Color', colors(1), 'LineStyle','-', 'LineWidth',2, 'DisplayName', 'Ishii');


plot(ubd_camber(:, 1), ubd_camber(:, 2), 'LineStyle','--', 'LineWidth',1.5, 'Color',colors(2), 'DisplayName','Camber line UBD-5459')
plot(ishii_camber(:, 1), ishii_camber(:, 2), 'LineStyle','--', 'LineWidth',1.5, 'Color',colors(1), 'DisplayName','Camber line Ishii')

plot(ubd_x, ubd_y, 'Color', colors(2), 'LineStyle','-', 'LineWidth',2, HandleVisibility='off');
plot(ishii_x, ishii_y, 'Color', colors(1), 'LineStyle','-', 'LineWidth',2, 'HandleVisibility', 'off');


axis equal;
grid on;
legend('Location','southeast', 'Color', 'w', 'EdgeColor', [0.1, 0.1, 0.1], 'TextColor','k', 'Interpreter','latex')
set(gca, 'Color', 'white', 'Xcolor', 'k', 'YColor', 'k', 'FontSize', fontsize_general,'GridColor', [0.05, 0.05, 0.05], 'XMinorGrid', 'off', 'YMinorGrid', 'off');


xlabel('$x/c$', 'Interpreter','latex', 'FontSize', 22);
ylabel('$y/c$', 'Interpreter','latex', 'FontSize', 22);

%ylim([-0.25, 0.25])
set(gcf, 'Color', 'w')


%% ------------------------ ISHII vs UBD --------------------------
%% Data loading
data_ishii_free = load('ishii_data_free_transition.txt');
data_ubd = load('UBD5494_Data_Clean.txt');

figure (2), clf;
hold on
ylim([-0.5, 1.5]);
xlim([-2, 10]);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


plot(data_ishii_free(:,1), data_ishii_free(:,2), 'LineWidth', 2, ...
        'DisplayName', ["Ishii"], ...
        'Color', colors(1), 'Marker','o');

plot(data_ubd(:,1), data_ubd(:,2), 'LineWidth', 2, ...
        'DisplayName', ["UBD-5459"], ...
        'Color', colors(2), 'Marker','o');


legend('Location', 'southeast', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_l$ [-]', 'Interpreter', 'latex', 'FontSize', 22);

grid on;
hold off
%%
figure (3), clf;
hold on;
xlim([0, 0.2]);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');

plot(data_ishii_free(:,3), data_ishii_free(:,2), 'LineWidth', 2, ...
        'DisplayName', ["Ishii"], ...
        'Color', colors(1), 'Marker','o');

plot(data_ubd(:,3), data_ubd(:,2), 'LineWidth', 2, ...
        'DisplayName', ["UBD-5459"], ...
        'Color', colors(2), 'Marker','o');

legend('Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('$C_d$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_l$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
grid on;
hold off


%% --------------------------Xtr comparison ----------------------
%% Data loading
data_ishii_free = load('ishii_data_free_transition.txt');
data_ishii_0_1 = load('ishii_data_0.1_transition.txt');
data_ishii_0_5 = load('ishii_data_0.5_transition.txt');

figure (4), clf;
hold on
ylim([-0.2, 1.2]);
yticks(-0.2:0.2:1.2)
xlim([-2, 10]);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


plot(data_ishii_free(:,1), data_ishii_free(:,2), 'LineWidth', 2, ...
        'DisplayName', ["No turbulator"], ...
        'Color', colors(1), 'Marker','o');
hold on
plot(data_ishii_0_1(:,1), data_ishii_0_1(:,2), 'LineWidth', 2, ...
        'DisplayName', ["$Xtr_{max}=0.1$"], ...
        'Color', colors(2), 'Marker','o');
plot(data_ishii_0_5(:,1), data_ishii_0_5(:,2), 'LineWidth', 2, ...
        'DisplayName', ["$Xtr_{max}=0.5$"], ...
        'Color', colors(3),'Marker','o');
legend('Location', 'southeast', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_l$ [-]', 'Interpreter', 'latex', 'FontSize', 22);

grid on;
hold off


%%
figure (5), clf;
hold on

xlim([0, 0.16]);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');

plot(data_ishii_free(:,3), data_ishii_free(:,2), 'LineWidth', 2, ...
        'DisplayName', ["No turbulator"], ...
        'Color', colors(1), 'Marker','o');
plot(data_ishii_0_1(:,3), data_ishii_0_1(:,2), 'LineWidth', 2, ...
        'DisplayName', ["$Xtr_{max}=0.1$"], ...
        'Color', colors(2), 'Marker','o');
plot(data_ishii_0_5(:,3), data_ishii_0_5(:,2), 'LineWidth', 2, ...
        'DisplayName', ["$Xtr_{max}=0.5$"], ...
        'Color', colors(3),'Marker','o');
legend('Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('$C_d$ [-]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_l$ [-]', 'Interpreter', 'latex', 'FontSize', 22);


set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

grid on;
hold off

%% Finding 
alpha_lo_free = interp1(data_ishii_free(:,2), data_ishii_free(:,1), 0);
alpha_lo_0_1  = interp1(data_ishii_0_1(:,2), data_ishii_0_1(:,1), 0);
alpha_lo_0_5  = interp1(data_ishii_0_5(:,2), data_ishii_0_5(:,1), 0);
fprintf("Ishii free transition: %f\n", alpha_lo_free)
fprintf("Ishii 0.1 fixed trans: %f\n", alpha_lo_0_1 )
fprintf("Ishii 0.5 fixed trans: %f\n", alpha_lo_0_5 )




