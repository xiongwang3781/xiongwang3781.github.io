% -----------------------------------------------------------------------------
% Author: Xiong Wang
% Date:   2025-06-09
% Purpose: Generate deformation mechanism map for dry olivine at T=1223 K, P=1.5 GPa,
%          overlay experimental grain‐size vs. stress data from multiple locations.
% -----------------------------------------------------------------------------

% Clear workspace and figures
clearvars;
close all;

%% Constants
R       = 8.314;        % Universal gas constant [J/(mol·K)]
T       = 1223;         % Temperature [K]
P       = 1.5;          % Pressure [GPa] (used in activation energy)

% Flow‐law parameters (Warren & Hansen, 2023)
% Dislocation creep
A_dis   = 1.1e5;        
n_dis   = 3.5;          
Q_dis   = (530 + 15*P)*1e3;  

% Diffusion creep
A_diff  = 10^7.6;       
n_diff  = 1;            
p_diff  = 3;            
Q_diff  = (375 + 1.5*6)*1e3;  

% Grain-boundary sliding (GBS)
A_gbs   = 10^4.8;       
n_gbs   = 2.9;          
p_gbs   = 0.7;          
Q_gbs   = (445 + 15*P)*1e3;  

% Generate grid of grain sizes and stresses
grain_size = logspace(0, 4, 400);    % Grain size [µm]
stress     = logspace(-1, 2, 400);   % Differential stress [MPa]
[GS, SS]   = meshgrid(grain_size, stress);

% Compute strain rates for each mechanism
rate_dis  = A_dis  * SS.^n_dis  .* exp(-Q_dis./(R.*T));
rate_diff = A_diff * SS.^n_diff .* GS.^(-p_diff) .* exp(-Q_diff./(R.*T));
rate_gbs  = A_gbs  * SS.^n_gbs  .* GS.^(-p_gbs) .* exp(-Q_gbs./(R.*T));

% Identify dominant mechanism at each grid point
% 1 = Dislocation, 2 = Diffusion, 3 = GBS
dominant = zeros(size(rate_dis));
dominant(rate_dis  >= max(rate_diff, rate_gbs)) = 1;
dominant(rate_diff >= max(rate_dis, rate_gbs)) = 2;
dominant(rate_gbs  >= max(rate_dis, rate_diff)) = 3;

% Mask out non-dominant rates for contour overlays
rate_dis(dominant~=1)  = NaN;
rate_diff(dominant~=2) = NaN;
rate_gbs(dominant~=3)  = NaN;

% Plot deformation mechanism map
figure('Color','w');
% Filled contours for dominant mechanisms
contourf(GS, SS, dominant, [0.5 1.5 2.5 3.5], 'LineColor','k','LineWidth',1.5);
colormap([0.4 0.76 0.65; 0.99 0.55 0.39; 0.55 0.63 0.79]);
hold on;

% Strain‐rate contour levels
levels = logspace(-18, -9, 10);

% Overlay strain‐rate contours (dashed)
contour(GS, SS, rate_dis,  levels, 'LineColor','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off');
[Cdiff, hDiff] = contour(GS, SS, rate_diff, levels, 'LineColor','k', 'LineStyle','--', 'LineWidth',1);
contour(GS, SS, rate_gbs,  levels, 'LineColor','k', 'LineStyle','--', 'LineWidth',1, 'HandleVisibility','off');

% Label diffusion‐creep contours only
clabel(Cdiff, hDiff, 'FontSize',12, 'LabelSpacing',1000);

% Annotate regions
text(4e3, 55, 'Dislocation', 'FontSize',12, 'HorizontalAlignment','center');
text(5,    1, 'Diffusion',   'FontSize',12, 'HorizontalAlignment','center');
text(5e2,  20, 'GBS',         'FontSize',12, 'HorizontalAlignment','center');

% Log scales and labels
set(gca, 'XScale','log', 'YScale','log');
xlim([1, 1e4]); ylim([0.1, 1e2]);
xlabel('Grain Size (\mum)',            'FontSize',14);
ylabel('Differential Stress (MPa)',    'FontSize',14);
title(sprintf('Dry Olivine, T = %d K, P = %.1f GPa', T, P), 'FontSize',16);

%% Load and overlay experimental data
data = readtable('Grainsize-Stress.txt', 'Delimiter','\t');
locations = unique(data.Location);
markers   = {'o','^','s'};
markerSize = 80;

hold on;
for i = 1:numel(locations)
    loc = locations{i};
    subset = data(strcmp(data.Location, loc), :);
    scatter(subset.GrainSize, subset.Stress, markerSize, ...
            'Marker', markers{i}, 'MarkerEdgeColor','r', ...
            'DisplayName', loc);
end
hold off;

set(gca, 'Layer','top');

% Save figure
print('Dry_olivine_deformation_map.png','-dpng','-r300');
