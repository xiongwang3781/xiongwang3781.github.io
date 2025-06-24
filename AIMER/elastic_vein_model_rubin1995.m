% Author: Xiong Wang
% Email: xiong.wang3781@outlook.com
% Date: 2025/06/24
%%
% === 参数设置 ===
E = 80;      % Young's modulus (GPa), typical crustal rock
nu = 0.25;   % Poisson's ratio
M = E / (1 - nu^2);  % Plane strain modulus (GPa)

% 扫描 w/(2L) 区间
ratio = linspace(0.005, 0.02, 100);  % w/(2L)，无量纲
dP = M .* ratio;                     % 差应力 (GPa)

% 美化参数全局设置
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultLegendFontName', 'Times New Roman');
set(0, 'DefaultLegendFontSize', 12);

% === 绘图 ===
fig = figure('Color','w', 'Position', [100, 100, 600, 400]);
plot(ratio, dP, 'r-', 'LineWidth', 2, 'DisplayName', 'E=80 GPa, \nu=0.25');
hold on;

% 示例点：w=1cm, L=38cm
w_cm = 1.0; L_cm = 38;
example_ratio = (w_cm/2)/L_cm;
example_dP = M * example_ratio;
plot(example_ratio, example_dP, 'ko', 'MarkerSize',10, 'MarkerFaceColor','k', ...
    'DisplayName','Field Example');

% 标签和标题
xlabel('w/(2L)');
ylabel('P_{melt} - |\sigma_x| (GPa)');

title('$P_{melt} - |\sigma_x| = M \frac{w}{2L} (Rubin 1995)$', ...
    'Interpreter', 'latex');


% 坐标轴及边框美化
set(gca, 'LineWidth', 1.2, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'Layer', 'top');
set(gca, 'TickDir', 'in', 'TickLength', [0.015 0.015]);

% legend
legend('show', 'Location', 'NorthWest', 'EdgeColor', [0 0 0], 'Box', 'off');

% 高分辨率保存
print(fig, 'differential_pressure_vs_ratio', '-dpng', '-r600');

% exportgraphics(fig, 'differential_pressure_vs_ratio.png', 'Resolution', 600);
% exportgraphics(fig, 'differential_pressure_vs_ratio.pdf', 'ContentType','vector');
