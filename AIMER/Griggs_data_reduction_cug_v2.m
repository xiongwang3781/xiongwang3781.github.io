%% Script for Griggs Raw Data Reduction
% Author: Xiong Wang
% Date: 2024/05/17
% Function: Find hitpoint, correct strain and stress, calculate flow stress
% Optimized by: Xiong Wang
% Date: 2025/04/01
clc;
clear;
close all;

%% Define Constants
stiffness_WC = 19; %um/kN
sample_initial_diameter = 3000;  % um, initial sample diameter
piston_diameter = 4763;          % um, piston diameter
piston_area = pi * (piston_diameter / 2)^2; %;
%% Read Excel File
file_path = 'C:\Users\xiong\OneDrive\Documents\Paper9\Griggs_data\GB013.xlsx';
[file_dir, file_name, ~] = fileparts(file_path); % Extract file directory and name
[data, txt, ~] = xlsread(file_path);                % Read Excel data
column_names = txt(3, :);           % 所有列名
% Extract experimental data
sample_initial_length = data(1, 2);      % um, initial sample length, 把样品初始长度数据放置在excel的这个单元格
time = data(5:end, 2);                    % NUM data
LVDT1 = data(5:end, 3);                  % LVDT1 data (um)
differential_stress = data(5:end, 7);    % Differential stress data (MPa)
pressure = data(5:end, 6); 
load = differential_stress * piston_area * 1e-9; % 单位：kN
temp_wc = data(5:end, 22);
temp_tc1 = data(5:end, 23);
temp_tc2 = data(5:end, 24);

%% Filter Data Where Differential Stress > 0
if ismember(file_name, {'GA422', 'GA429'})
    valid_idx = (differential_stress > 0); 
else
    valid_idx = (differential_stress > 0) & (time > 30000); 
end

time_valid = time(valid_idx);
LVDT1_valid = LVDT1(valid_idx);
differential_stress_valid = differential_stress(valid_idx);
pressure_valid = pressure(valid_idx);
load_valid = load(valid_idx);
temp_wc_valid = temp_wc(valid_idx);
temp_tc1_valid = temp_tc1(valid_idx);
temp_tc2_valid = temp_tc2(valid_idx);
%% Plot NUM vs Differential Stress and Select Range
figure('Units', 'inches', 'Position', [1, 1, 12, 6], 'Color', 'white');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
plot(time_valid, differential_stress_valid, 'o-', 'MarkerSize', 3, 'LineWidth', 1.0, 'Color','k','MarkerEdgeColor','k','MarkerFaceColor','k');
%plot(NUM_valid, differential_stress_valid, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5);
title('Please click two points on the graph to select the range', 'FontSize', 14);
xlabel('Time (s)');
ylabel('Differential Stress (MPa)');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');

disp('Please select two points on the NUM vs Differential Stress plot to define the range.');
[x, ~] = ginput(2); % User selects two points

% Find closest actual data points
[~, idx1] = min(abs(time_valid - x(1)));
[~, idx2] = min(abs(time_valid - x(2)));
xmin = min(time_valid(idx1), time_valid(idx2));
xmax = max(time_valid(idx1), time_valid(idx2));
disp(['Selected range: ', num2str(xmin), ' to ', num2str(xmax)]);

%% Filter Data Within the Selected Range
filtered_idx = (time_valid >= xmin) & (time_valid <= xmax);
filtered_time = time_valid(filtered_idx);
filtered_LVDT1 = LVDT1_valid(filtered_idx);
filtered_stress = differential_stress_valid(filtered_idx);
filtered_pressure = pressure_valid(filtered_idx);
filtered_load = load_valid(filtered_idx);
filtered_temp_wc = temp_wc_valid(filtered_idx);
filtered_temp_tc1 = temp_tc1_valid(filtered_idx);
filtered_temp_tc2 = temp_tc2_valid(filtered_idx);

if isempty(filtered_time)
    error('No data points in the selected range.');
end

%% Plot LVDT1 vs Differential Stress and Select Two Ranges for Linear Fits
figure('Units', 'inches', 'Position', [1, 1, 5, 4], 'Color', 'white');
plot(filtered_LVDT1, filtered_stress, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5, 'Color','k','MarkerEdgeColor','k','MarkerFaceColor','k');
%plot(filtered_LVDT1, filtered_stress, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5);
title('Please click four points on the graph', 'FontSize', 14);
xlabel('LVDT1 (um)');
ylabel('Differential Stress (MPa)');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k','XDir','reverse');

disp('Please select four points on the LVDT1 vs Differential Stress plot to define two ranges for linear fits.');
[x, ~] = ginput(4); % User selects four points

% Find closest actual data points
[~, idx1] = min(abs(filtered_LVDT1 - x(1)));
[~, idx2] = min(abs(filtered_LVDT1 - x(2)));
[~, idx3] = min(abs(filtered_LVDT1 - x(3)));
[~, idx4] = min(abs(filtered_LVDT1 - x(4)));
selected_xmin1 = min(filtered_LVDT1(idx1), filtered_LVDT1(idx2));
selected_xmax1 = max(filtered_LVDT1(idx1), filtered_LVDT1(idx2));
selected_xmin2 = min(filtered_LVDT1(idx3), filtered_LVDT1(idx4));
selected_xmax2 = max(filtered_LVDT1(idx3), filtered_LVDT1(idx4));
disp(['Range 1: ', num2str(selected_xmin1), ' to ', num2str(selected_xmax1)]);
disp(['Range 2: ', num2str(selected_xmin2), ' to ', num2str(selected_xmax2)]);

%% Filter and Fit the Two Selected Ranges
range1_idx = (filtered_LVDT1 >= selected_xmin1) & (filtered_LVDT1 <= selected_xmax1);
range2_idx = (filtered_LVDT1 >= selected_xmin2) & (filtered_LVDT1 <= selected_xmax2);

if isempty(filtered_LVDT1(range1_idx)) || isempty(filtered_LVDT1(range2_idx))
    error('No data points in one or both selected ranges.');
end

range1_LVDT1 = filtered_LVDT1(range1_idx);
range1_stress = filtered_stress(range1_idx);
range2_LVDT1 = filtered_LVDT1(range2_idx);
range2_stress = filtered_stress(range2_idx);

% Linear regression
p1 = polyfit(range1_LVDT1, range1_stress, 1);
p2 = polyfit(range2_LVDT1, range2_stress, 1);

% Calculate intersection point
x_intersect = (p2(2) - p1(2)) / (p1(1) - p2(1));
y_intersect = polyval(p1, x_intersect);

%% Plot LVDT1 vs Differential Stress with Linear Fits
figure('Units', 'inches', 'Position', [1, 1, 5, 4], 'Color', 'white');
plot(filtered_LVDT1, filtered_stress, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5, 'Color','k','MarkerEdgeColor','k','MarkerFaceColor','k'); hold on;

%plot(filtered_LVDT1, filtered_stress, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5); hold on;
plot(range1_LVDT1, polyval(p1, range1_LVDT1), 'b-', 'LineWidth', 2);
plot(range2_LVDT1, polyval(p2, range2_LVDT1), 'r-', 'LineWidth', 2);

% Extend fit lines
x_fit_1_ext = linspace(min(range1_LVDT1) - 400, max(range1_LVDT1), 100);
x_fit_2_ext = linspace(min(range2_LVDT1), max(range2_LVDT1) + 200, 100);
plot(x_fit_1_ext, polyval(p1, x_fit_1_ext), 'b--', 'LineWidth', 1);
plot(x_fit_2_ext, polyval(p2, x_fit_2_ext), 'r--', 'LineWidth', 1);

% Annotate intersection point
text(x_intersect + 10, y_intersect, sprintf('Hitpoint: (%0.1f, %0.1f)', x_intersect, y_intersect), ...
    'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%title('LVDT1 vs Differential Stress with Linear Fits', 'FontSize', 14);
xlabel('LVDT1 (um)');
ylabel('Differential Stress (MPa)');

set(gca,'XDir','reverse','Box','on','LineWidth',1.2,'XColor','k','YColor','k');
legend('Raw data', 'Fit 1', 'Fit 2', 'Location', 'southeast');

% Save figure
output_figure_path = fullfile(file_dir, [file_name, '_hitpoint_plot.png']);
print(gcf, output_figure_path, '-dpng', '-r300');
hold off;
disp(['Intersection point: (', num2str(x_intersect), ', ', num2str(y_intersect), ')']);

%% Calculate Strain and Stress After Intersection Point

post_intersection_idx = filtered_LVDT1 <= x_intersect;
post_intersection_load = filtered_load(post_intersection_idx);
post_intersection_LVDT1 = filtered_LVDT1(post_intersection_idx);
post_intersection_time = filtered_time(post_intersection_idx);
post_intersection_stress = filtered_stress(post_intersection_idx);
load_relative = post_intersection_load - post_intersection_load(1);
displacement = x_intersect - post_intersection_LVDT1;
displacement_corrected = x_intersect - post_intersection_LVDT1 - load_relative.*stiffness_WC; % corrected for WC stiffness
strain = -log(1 - displacement_corrected / sample_initial_length);
friction = y_intersect;
stress = (1 - displacement_corrected / sample_initial_length) .* ...
         (post_intersection_stress - friction) * ...
         (piston_diameter / sample_initial_diameter)^2;

% Smooth stress data
window_size = 40; % Adjust based on experiment speed, e.g., 400 for slow experiments
stress_smoothed = movmean(stress, window_size);

%% plot temp and pressure
figure('Units','inches','Position',[1,1,5,4],'Color','white');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
post_intersection_pressure = filtered_pressure(post_intersection_idx);
post_intersection_temp_wc = filtered_temp_wc(post_intersection_idx);
post_intersection_temp_tc1 = filtered_temp_tc1(post_intersection_idx);
post_intersection_temp_tc2 = filtered_temp_tc2(post_intersection_idx);
time_offset = post_intersection_time - post_intersection_time(1);

plot(time_offset, post_intersection_temp_wc, 'r-', 'LineWidth', 4.0);
hold on;
plot(time_offset, post_intersection_temp_tc1, 'b-', 'LineWidth', 1.5);
plot(time_offset, post_intersection_temp_tc2, 'ko', 'MarkerSize', 1, 'LineWidth', 1.5);
%plot(time_offset, post_intersection_temp_tc2, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature (°C)');
ylim([850 1050]);
yyaxis right
hold on;
plot(time_offset, post_intersection_pressure, 'k-', 'LineWidth', 1.5);
xlim([0, time_offset(end)]);
ylim([800 1500]);
yticks(500:100:1500);     
ylabel('Pressure (MPa)');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');

legend('Working TC','TC1','TC2','Pressure', 'Location','southeast');
output_figure_path = fullfile(file_dir, [file_name, '_temp_pressure.png']);
print(gcf, output_figure_path, '-dpng', '-r300'); 
%% Plot Strain vs Smoothed Stress
fig = figure('Units', 'inches', 'Position', [1, 1, 5, 4], 'Color', 'white');
plot(strain, stress_smoothed, 'k-', 'LineWidth', 1.5);

% Annotate peak stress
[max_stress_value, max_stress_index] = max(stress_smoothed);
text(strain(max_stress_index) + 0.001, max_stress_value, ...
    sprintf('Peak Stress: %0.2f MPa', max_stress_value), ...
    'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
%title('Stress-Strain Curve', 'FontSize', 14);
xlabel('Strain');
ylabel('Stress (MPa)');
%grid on;
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');
% Save figure
output_figure_path = fullfile(file_dir, [file_name, '_stress_strain_plot.png']);
print(gcf, output_figure_path, '-dpng', '-r300');

%% Save Strain and Smoothed Stress Data
output_data_path = fullfile(file_dir, [file_name, '_stress_strain_data.txt']);
fid = fopen(output_data_path, 'w');
fprintf(fid, 'Strain\tStress_smoothed\n');
for i = 1:length(strain)
    fprintf(fid, '%f\t%f\n', strain(i), stress_smoothed(i));
end
fclose(fid);
disp(['Data saved to: ', output_data_path]);

%%
figure('Units','inches','Position',[1,1,15,4],'Color','white');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
% Shift time so it starts at 0
time_offset = post_intersection_time - post_intersection_time(1);
window_size = 10;  % can adjust based on preference or sampling rate
displacement_smoothed = movmean(displacement, window_size, 'Endpoints', 'fill');
% Left subplot (1): Time vs Displacement
subplot(1,3,1);
plot(time_offset, displacement_smoothed, 'k-','LineWidth',1.5);
xlim([0, time_offset(end)]);
xlabel('Time (s)');
ylabel('Piston disp. (\mum)');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');

subplot(1,3,2);
% Calculate displacement rate (velocity)
velocity = diff(displacement_smoothed) ./ diff(time_offset);  % displacement rate
% Time midpoints for velocity
time_mid = 0.5*(time_offset(1:end-1) + time_offset(2:end));
% Middle subplot (2): Time vs Displacement Rate
plot(time_mid, velocity, 'k-','LineWidth',1.5);

% Smooth velocity (movmean)
window_size = 50;  % can adjust based on preference or sampling rate
velocity_smoothed = movmean(velocity, window_size);
hold on;
plot(time_mid, velocity_smoothed, 'r-','LineWidth',1.5);
xlim([0, time_mid(end)]);
xlabel('Time (s)');
ylabel('Piston vel. (\mum/s)');
legend('Raw velocity','Smoothed velocity','Location','northeast');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');

% Calculate strain rate = displacement rate / sample_initial_length
% sample_initial_length is in microns, so strain_rate is in 1/s if time is in s
strain_rate = velocity / sample_initial_length;

% Smooth strain rate
strain_rate_smoothed = movmean(strain_rate, window_size);

% Right subplot (3): Time vs Strain Rate
subplot(1,3,3);
plot(time_mid, strain_rate, 'k-','LineWidth',1.5);
hold on;
plot(time_mid, strain_rate_smoothed, 'r-','LineWidth',1.5);
xlim([0, time_mid(end)]);
xlabel('Time (s)');
ylabel('Strain rate (1/s)');
legend('Raw strain rate','Smoothed strain rate','Location','northeast');
set(gca,'Box','on','LineWidth',1.2,'XColor','k','YColor','k');

% Save figure
output_figure_path = fullfile(file_dir, [file_name, '_time_displacement_velocity_strainrate.png']);
print(gcf, output_figure_path, '-dpng', '-r300');
disp(['Time-displacement, velocity, and strain rate figure saved to: ', output_figure_path]);

