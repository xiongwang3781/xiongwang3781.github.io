% Author: Xiong Wang
% Email: xiong.wang3781@outlook.com
% Date: 2025/06/24

%%
clear; clc;
% ==== PLOT RESULT (YEARS) ====
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 16);
set(0, 'DefaultLegendFontName', 'Times New Roman');
set(0, 'DefaultLegendFontSize', 14);
%%

clear; clc;

% =======================
% PARAMETERS
% =======================
D = 2e-17;      % Diffusion coefficient of Ca in olivine (cm^2/s)
R = 100e-4;     % Crystal radius (100 microns, converted to cm)
Nx = 101;       % Number of spatial grid points
dx = R/(Nx-1);  % Spatial step (cm)
x = linspace(0, R, Nx);    % Radial position array (cm)
x_um = x * 1e4;            % Radial position array (microns)

% =======================
% 1. GENERATE SYNTHETIC 'EXPERIMENTAL PROFILE'
%    (Sigmoid shape, core-to-rim decrease, simulates a real EPMA traverse)
% =======================
center = 40;    % Center position of the diffusion front (microns)
width = 5;      % Slope width (microns)
c_exp = 0.1 + 0.9 ./ (1 + exp((x_um-center)/width));  % Sigmoid profile: core-high, rim-low
c_exp = c_exp + 0.01*randn(size(c_exp));             % Optional: add small noise to mimic measurement

% =======================
% 2. INITIAL PROFILE FOR DIFFUSION MODEL
%    (Sharp step: Ca=1 for core, Ca=0 for rim, before diffusion)
% =======================
c0 = zeros(Nx,1);        % Start with zero everywhere
c0(x < center*1e-4) = 1; % For radii < center (in cm), Ca=1; elsewhere, Ca=0

% =======================
% 3. AUTOMATED FITTING LOOP
%    (Test a series of diffusion times, find the best match to 'c_exp')
% =======================
test_times = logspace(10, 11, 100); % Test 100 candidate diffusion times (s)
years = test_times / (365.25*24*3600); %
disp(['Min: ', num2str(min(years), '%.0f'), ' years, Max: ', num2str(max(years), '%.0f'), ' years'])

best_err = Inf;        % Initialize lowest error
best_profile = [];     % Placeholder for best-fit profile
best_time = 0;         % Placeholder for best-fit time

for t_idx = 1:length(test_times)
    t = test_times(t_idx);              % Candidate diffusion time (s)
    Nt = ceil(t / (0.5 * dx^2 / D));    % Number of time steps, ensures stability
    dt = t / Nt;                        % Time step size (s)
    c = c0;                             % Start with initial (step) profile
    
    % Finite difference solution (explicit method)
    for n = 1:Nt
        c_new = c;
        for i = 2:Nx-1
            c_new(i) = c(i) + D * dt / dx^2 * (c(i+1) - 2*c(i) + c(i-1));
        end
        % No-flux (Neumann) boundary condition at core and rim
        c_new(1) = c_new(2);
        c_new(end) = c_new(end-1);
        c = c_new;
    end
    
    model_profile = c';                             % Simulated profile at this diffusion time
    err = mean((model_profile - c_exp).^2);         % Mean squared error vs. 'experimental' profile
    
    % Keep best fit
    if err < best_err
        best_err = err;
        best_profile = model_profile;
        best_time = t;
    end
end

% =======================
% 4. RESULT DISPLAY (IN YEARS)
%    Plot measured vs. best-fit profiles and report best-fitting diffusion time
% =======================
seconds_in_year = 365.25 * 24 * 3600;           % Conversion factor: seconds to years
best_time_yr = best_time / seconds_in_year;     % Best-fit diffusion time (years)

figure('Color','w', 'Position', [100, 100, 600, 400]); 
hold on;
plot(x_um, c_exp, 'ko-', 'LineWidth',2, 'DisplayName','This study');
plot(x_um, best_profile, 'r-', 'LineWidth',2, 'DisplayName','Best fitting');
xlabel('Radius (\mum)');
ylabel('Ca content');
legend show; box on;
title(['Best fitting diffusion time t = ', num2str(best_time_yr, '%.0f'), ' years']);

set(gca, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
set(gca, 'TickDir', 'in', 'TickLength', [0.015 0.015]);

legend('show', 'Location', 'best', 'EdgeColor', [0 0 0], 'Box', 'off');

disp(['Best fitting diffusion time: ', num2str(best_time), ' s (', ...
    num2str(best_time/3600, '%.2f'), ' hours, ', ...
    num2str(best_time_yr, '%.4f'), ' years)']);

% ===== END =====
% This script generates a synthetic Ca diffusion profile (core-high to rim-low),
% fits the best diffusion time using a 1D finite-difference model,
% and reports the optimal duration in seconds, hours, and years.

%%
clear; clc;

% ==== PARAMETERS ====
D = 2e-17;      % Diffusion coefficient of Ca in olivine (cm^2/s)
R = 100e-4;     % Crystal radius (100 microns, converted to cm)
Nx = 101;       % Number of spatial grid points
dx = R/(Nx-1);  % Spatial step (cm)
x = linspace(0, R, Nx);    % Radial position array (cm)
x_um = x * 1e4;            % Radial position array (microns)

% ==== SYNTHETIC "EPMA" PROFILE (sigmoid, core low, rim high) ====
c_core = 0.02;     % Initial Ca in olivine core (e.g., wt%)
c_melt = 0.18;     % Ca at rim (melt boundary, e.g., wt%)
center = 80;       % Position of transition, microns
width = 8;         % Slope width (microns)
c_exp = c_core + (c_melt-c_core) ./ (1 + exp(-(x_um-center)/width)); % sigmoid
c_exp = c_exp + 0.002*randn(size(c_exp)); % small noise

% ==== INITIAL PROFILE: UNIFORM CORE Ca ====
c0 = c_core * ones(Nx,1);

% ==== AUTOMATIC DIFFUSION TIME FITTING ====
test_times = logspace(10, 11.7, 200); 
years = test_times / (365.25*24*3600);
disp(['Min: ', num2str(min(years), '%.0f'), ' years, Max: ', num2str(max(years), '%.0f'), ' years'])

best_err = Inf; best_profile = []; best_time = 0;
plot_profiles = []; plot_times = []; % for visualization

for t_idx = 1:length(test_times)
    t = test_times(t_idx);
    Nt = ceil(t / (0.4 * dx^2 / D));  % number of time steps for stability
    dt = t / Nt;
    c = c0;
    for n = 1:Nt
        c_new = c;
        for i = 2:Nx-1
            c_new(i) = c(i) + D * dt / dx^2 * (c(i+1) - 2*c(i) + c(i-1));
        end
        % Boundary conditions
        c_new(1) = c_new(2);      % no flux at core
        c_new(end) = c_melt;      % fixed Ca at rim
        c = c_new;
    end
    model_profile = c';
    err = mean((model_profile - c_exp).^2); % mean squared error
    if err < best_err
        best_err = err;
        best_profile = model_profile;
        best_time = t;
    end
    % Save some profiles for visualization
    if ismember(t_idx, [1 50 100 150 200])
        plot_profiles = [plot_profiles; model_profile];
        plot_times = [plot_times; t];
    end
end



seconds_in_year = 365.25 * 24 * 3600;
best_time_yr = best_time / seconds_in_year;

fig = figure('Color','w', 'Position', [100, 100, 600, 400]);
hold on;
plot(x_um, c_exp, 'ko-', 'LineWidth',2, 'MarkerFaceColor','w', 'DisplayName','Data');
plot(x_um, best_profile, 'r-', 'LineWidth',2, 'DisplayName','Best fit');

colors = lines(length(plot_times));
for i = 1:length(plot_times)
    plot(x_um, plot_profiles(i,:), '--', 'Color',colors(i,:), ...
        'LineWidth', 1.2, ...
        'DisplayName',['t=',num2str(plot_times(i)/seconds_in_year,'%.0f'),' yr']);
end

xlabel('Distance (\mum)');
ylabel('Ca content');
title(['Best fitting diffusion time t = ', num2str(best_time_yr, '%.0f'), ' years']);

set(gca, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
set(gca, 'TickDir', 'in', 'TickLength', [0.015 0.015]);

legend('show', 'Location', 'NorthWest', 'EdgeColor', [0 0 0], 'Box', 'off');

disp(['Best fitting diffusion time: ', num2str(best_time), ' s (', ...
    num2str(best_time_yr, '%.4f'), ' years)']);

% 保存为高分辨率PNG（600 dpi）：
print(fig, 'Ca_diffusion_profile', '-dpng', '-r600');

% % 保存为高分辨率TIFF（一般用在SCI期刊）：
% print(fig, 'Ca_diffusion_profile', '-dtiff', '-r600');
% 
% % 矢量PDF/EPS（适合所有印刷/排版）：
% print(fig, 'Ca_diffusion_profile', '-dpdf');   % PDF矢量
% print(fig, 'Ca_diffusion_profile', '-depsc');  % EPS矢量


