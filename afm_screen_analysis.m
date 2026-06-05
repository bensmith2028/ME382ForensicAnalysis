clear; close all; clc;

% initial scan parameters
DATA_DIR = 'RawData';
SCAN_UM  = 25; % scan size in micrometers
N        = 256; % pixels per side

% sample definitions
tags = { ...
    'anti1-lateral-mode_060126112724', ...
    'anti2-lateral-mode_060126115108', ...
    'non1-lateral-mode_060126120811',  ...
    'non2-lateral-mode_060126122209'   };

labels = {'Anti-Glare #1', 'Anti-Glare #2', ...
          'Tempered Glass #1', 'Tempered Glass #2'};

% colors, anti-glare: blue, tempered glass: red
clr = {'#1f77b4', '#aec7e8', '#d62728', '#f4a582'};

nSamp = numel(tags);
x_um  = linspace(0, SCAN_UM, N);
y_um  = linspace(0, SCAN_UM, N);

% load/orient data
% forward scan direction + height scan: right-to-left per row, bottom to top.
%
% backward scan direction: left-to-right per row, bottom to top.
%
% friction average: at each spatial location FRW ≈ +F, BKW ≈ −F, average
% eliminate topography crosstalk
HEIGHT = cell(1, nSamp);
FRW    = cell(1, nSamp);
BKW    = cell(1, nSamp);
FRIC   = cell(1, nSamp);

for i = 1:nSamp
    H_raw  = load_channel(DATA_DIR, tags{i}, 'SIG_HEIGHT_SENSOR_FRW', N);
    F_raw  = load_channel(DATA_DIR, tags{i}, 'SIG_L_R_FRW',           N);
    B_raw  = load_channel(DATA_DIR, tags{i}, 'SIG_L_R_BKW',           N);

    HEIGHT{i} = fliplr(H_raw); % right-to-left to left-to-right
    FRW{i}    = fliplr(F_raw);
    BKW{i}    = B_raw; % already left-to-right

    HEIGHT{i} = plane_level(HEIGHT{i}); % remove linear tilt

    FRIC{i} = (FRW{i} - BKW{i}) / 2;
end

% step-corrected height: detect inter-row jumps via MAD-based threshold
% and apply cumulative offsets so each row's midpoint aligns with the prior row.
HEIGHT_CORR = HEIGHT;
for i = 1:nSamp
    HEIGHT_CORR{i} = correct_row_jumps(HEIGHT{i});
end

% figure 1: topography: 2D height maps
figure('Name', 'Topography 2-D', 'Position', [30 30 1200 620]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nSamp
    nexttile;
    imagesc(x_um, y_um, HEIGHT{i}); axis image xy;
    colormap(gca, 'hot');
    cb = colorbar; cb.Label.String = 'Height (nm)';
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    title(labels{i}, 'FontSize', 11);
end
sgtitle('Surface Topography — Height Channel (2-D)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% figure 2: topography: 3D height maps
figure('Name', 'Topography 3-D', 'Position', [50 50 1200 620]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nSamp
    nexttile;
    surf(x_um, y_um, HEIGHT{i}, 'EdgeColor', 'none');
    colormap(gca, 'hot'); shading interp;
    lighting gouraud; camlight headlight;
    xlabel('x  (\mum)'); ylabel('y  (\mum)'); zlabel('Height (nm)');
    title(labels{i}, 'FontSize', 11);
    view(-40, 28);
end
sgtitle('Surface Topography — Height Channel (3-D)', ...
    'FontSize', 13, 'FontWeight', 'bold');

%  figure 3: topography: 3D height maps (step-corrected)
figure('Name', 'Topography 3-D (Step-Corrected)', 'Position', [70 70 1200 620]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nSamp
    nexttile;
    surf(x_um, y_um, HEIGHT_CORR{i}, 'EdgeColor', 'none');
    colormap(gca, 'hot'); shading interp;
    lighting gouraud; camlight headlight;
    xlabel('x  (\mum)'); ylabel('y  (\mum)'); zlabel('Height (nm)');
    title(labels{i}, 'FontSize', 11);
    view(-40, 28);
end
sgtitle('Surface Topography — Step-Corrected (3-D)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% figure 4: friction maps
figure('Name', 'Friction Maps', 'Position', [80 80 1400 820]);
t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

row_data   = {FRW, BKW, FRIC};
row_titles = {'Friction — Forward (V)', ...
              'Friction — Backward (V)', ...
              'Friction — Average (V)'};

for row = 1:3
    for i = 1:nSamp
        nexttile;
        imagesc(x_um, y_um, row_data{row}{i}); axis image xy;
        if row < 3
            % diverging colormap centred at zero for signed channels
            colormap(gca, bwr_cmap(256));
            v = max(abs(row_data{row}{i}(:))) * 1.05;
            clim([-v, v]);
        else
            colormap(gca, 'parula');
        end
        cb = colorbar; cb.Label.String = 'Signal (V)';
        xlabel('x  (\mum)'); ylabel('y  (\mum)');
        title(labels{i}, 'FontSize', 10);
    end
end

for row = 1:3
    nexttile((row-1)*4 + 1);
    ylabel(sprintf('%s\ny (\\mum)', row_titles{row}), 'FontSize', 10);
end
sgtitle('Friction Channel Maps', 'FontSize', 13, 'FontWeight', 'bold');

% figure 5: cross-section profiles (horizontal midline
figure('Name', 'Height Profiles', 'Position', [100 100 1200 650]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
mid = round(N / 2);
for i = 1:nSamp
    nexttile;
    plot(x_um, HEIGHT_CORR{i}(mid, :), 'Color', clr{i}, 'LineWidth', 1.2);
    xlabel('x  (\mum)'); ylabel('Height (nm)');
    title(sprintf('%s  —  mid-row profile', labels{i}), 'FontSize', 11);
    grid on; xlim([0, SCAN_UM]);
end
sgtitle('Height Cross-Section Profiles (centre row)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% roughness metrics
fprintf('\n=== Roughness Metrics (Height Channel) ===\n');
fprintf('%-24s  %8s  %8s  %10s  %8s  %8s\n', ...
    'Sample', 'Rq (nm)', 'Ra (nm)', 'Rz (nm)', 'Rsk', 'Rku');
fprintf('%s\n', repmat('-', 1, 72));

Rq_vals   = zeros(1, nSamp);
Ra_vals   = zeros(1, nSamp);
Rz_vals   = zeros(1, nSamp);
Fric_mean = zeros(1, nSamp);
Fric_rms  = zeros(1, nSamp);
Fric_std  = zeros(1, nSamp);
Fric_cv   = zeros(1, nSamp);

for i = 1:nSamp
    h  = HEIGHT_CORR{i}(:);
    h  = h - mean(h);
    Rq = sqrt(mean(h.^2));
    Ra = mean(abs(h));
    Rz = max(h) - min(h);
    Rsk = mean(h.^3) / Rq^3;
    Rku = mean(h.^4) / Rq^4;
    Rq_vals(i) = Rq;
    Ra_vals(i) = Ra;
    Rz_vals(i) = Rz;
    fprintf('%-24s  %8.3f  %8.3f  %10.3f  %8.3f  %8.3f\n', ...
        labels{i}, Rq, Ra, Rz, Rsk, Rku);

    f = FRIC{i}(:);
    Fric_mean(i) = mean(f);
    Fric_rms(i)  = sqrt(mean(f.^2));
    Fric_std(i)  = std(f);
    Fric_cv(i)   = Fric_std(i) / abs(Fric_mean(i));
end
fprintf('\nAnti-Glare  mean Rq = %.3f nm  (avg of #1, #2)\n', ...
    mean(Rq_vals(1:2)));
fprintf('Temp. Glass mean Rq = %.3f nm  (avg of #1, #2)\n\n', ...
    mean(Rq_vals(3:4)));

fprintf('=== Friction Metrics (Average Channel) ===\n');
fprintf('%-24s  %12s  %12s\n', 'Sample', 'Mean (V)', 'RMS (V)');
fprintf('%s\n', repmat('-', 1, 52));
for i = 1:nSamp
    fprintf('%-24s  %12.4f  %12.4f\n', labels{i}, Fric_mean(i), Fric_rms(i));
end
fprintf('\nAnti-Glare  mean friction RMS = %.4f V\n', mean(Fric_rms(1:2)));
fprintf('Temp. Glass mean friction RMS = %.4f V\n\n', mean(Fric_rms(3:4)));

fprintf('=== Friction Variance Analysis ===\n');
fprintf('%-24s  %10s  %10s  %10s\n', 'Sample', 'Std (V)', 'Variance', 'CoV');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:nSamp
    fprintf('%-24s  %10.4f  %10.6f  %10.3f\n', ...
        labels{i}, Fric_std(i), Fric_std(i)^2, Fric_cv(i));
end
fprintf('\nAnti-Glare  mean friction Std = %.4f V  (variance = %.6f V^2)\n', ...
    mean(Fric_std(1:2)), mean(Fric_std(1:2).^2));
fprintf('Temp. Glass mean friction Std = %.4f V  (variance = %.6f V^2)\n', ...
    mean(Fric_std(3:4)), mean(Fric_std(3:4).^2));
fprintf('Variance ratio (Anti/Non) = %.2fx\n\n', ...
    mean(Fric_std(1:2).^2) / mean(Fric_std(3:4).^2));

% figure 6: structure function
px_nm  = SCAN_UM * 1e3 / N;
tau_nm = (1 : N-1) * px_nm;

figure('Name', 'Structure Function', 'Position', [120 120 1200 540]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

group_sf     = {1:2, 3:4};
group_sf_lbl = {'Anti-Glare', 'Tempered Glass'};

for g = 1:2
    nexttile; hold on;
    for i = group_sf{g}
        S = struct_func(HEIGHT_CORR{i});
        loglog(tau_nm, S, 'Color', clr{i}, 'LineWidth', 1.8, ...
            'DisplayName', labels{i});
    end
    xlabel('\tau  (nm)', 'FontSize', 12);
    ylabel("S'(\tau)  (nm^2)", 'FontSize', 12);
    title(group_sf_lbl{g}, 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 10);
    grid on; box on; hold off;
end
sgtitle("Structure Function  S'(\tau)  —  Height Channel", ...
    'FontSize', 13, 'FontWeight', 'bold');

% figure 7: roughness/friction bar charts
figure('Name', 'Roughness & Friction Comparison', 'Position', [150 150 1000 450]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

bar_x   = categorical(labels, labels);   % preserve order
rgb_mat = cell2mat(cellfun(@hex2rgb, clr', 'UniformOutput', false));

nexttile;
b1 = bar(bar_x, Rq_vals);
b1.FaceColor = 'flat';
b1.CData = rgb_mat;
ylabel('R_q  (nm)', 'FontSize', 12);
title('RMS Roughness  R_q', 'FontSize', 12);
grid on; box on;

nexttile;
b2 = bar(bar_x, Fric_rms);
b2.FaceColor = 'flat';
b2.CData = rgb_mat;
ylabel('Friction RMS  (V)', 'FontSize', 12);
title('Friction RMS (Average Channel)', 'FontSize', 12);
grid on; box on;

sgtitle('Material Comparison', 'FontSize', 13, 'FontWeight', 'bold');

% figure 8: local friction standard deviation
WIN = 16; % each pixel is standard deviation within WINxWIN window

figure('Name', 'Local Friction Std Dev', 'Position', [170 170 1200 620]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nSamp
    nexttile;
    Lsd = local_std_map(FRIC{i}, WIN);
    x_lsd = linspace(0, SCAN_UM, size(Lsd, 2));
    y_lsd = linspace(0, SCAN_UM, size(Lsd, 1));
    imagesc(x_lsd, y_lsd, Lsd); axis image xy;
    colormap(gca, 'hot');
    cb = colorbar; cb.Label.String = '\sigma_{friction}  (V)';
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    title([labels{i} '  —  local \sigma_{f}'], 'FontSize', 11);
end
sgtitle(sprintf('Local Friction Std Dev  (%d\\times%d px window)', WIN, WIN), ...
    'FontSize', 13, 'FontWeight', 'bold');

% figure 9: friction histograms
figure('Name', 'Friction Distributions', 'Position', [190 190 900 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

group_idx  = {1:2, 3:4};
group_name = {'Anti-Glare', 'Tempered Glass'};

for g = 1:2
    nexttile; hold on;
    for i = group_idx{g}
        f = FRIC{i}(:);
        histogram(f, 80, 'Normalization', 'probability', ...
            'FaceColor', clr{i}, 'FaceAlpha', 0.55, 'EdgeColor', 'none', ...
            'DisplayName', labels{i});
        xline(mean(f),   '--', 'Color', clr{i}, 'LineWidth', 1.4, ...
            'HandleVisibility', 'off');
        xline(mean(f)+std(f), ':', 'Color', clr{i}, 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
        xline(mean(f)-std(f), ':', 'Color', clr{i}, 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
    xlabel('Friction  (V)', 'FontSize', 12);
    ylabel('Probability', 'FontSize', 12);
    title(group_name{g}, 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 10);
    grid on; box on; hold off;
end
sgtitle('Friction Value Distributions  (dashed = mean, dotted = ±1\sigma)', ...
    'FontSize', 13, 'FontWeight', 'bold');

% local functions
% load data
function Z = load_channel(data_dir, tag, channel, N)
    fname = fullfile(data_dir, sprintf('%s.%s.FLT.txt', tag, channel));
    fid = fopen(fname, 'r');
    if fid == -1, error('Cannot open: %s', fname); end
    fgetl(fid);
    raw = fscanf(fid, '%f');
    fclose(fid);
    if numel(raw) ~= N * N
        error('Expected %d values in %s, got %d', N*N, fname, numel(raw));
    end
    Z = reshape(raw, N, N)';
end

% least-squares plane subtraction, remove sample tilt
function Z = plane_level(Z)
    [rows, cols] = size(Z);
    [X, Y] = meshgrid(1:cols, 1:rows);
    A = [X(:), Y(:), ones(rows * cols, 1)];
    coeffs = A \ Z(:);
    Z = Z - reshape(A * coeffs, rows, cols);
end

% S'(tau) = <[Z(x+tau) - Z(x)]^2>, averaged over all rows and x
function S = struct_func(Z)
    ncols = size(Z, 2);
    S = zeros(1, ncols - 1);
    for k = 1 : ncols - 1
        d = Z(:, k+1:end) - Z(:, 1:end-k);
        S(k) = mean(d(:).^2);
    end
end

% detect cross-row discontinuities and apply offset based on row midpoint
function Z_corr = correct_row_jumps(Z)
    nrows   = size(Z, 1);
    mid_col = round(size(Z, 2) / 2);
    half    = 5;   % average ±5 columns around centre for noise robustness
    cols    = max(1, mid_col-half) : min(size(Z,2), mid_col+half);
    mid_vals = mean(Z(:, cols), 2);   % representative midpoint value per row

    diffs     = diff(mid_vals);
    mad_val   = median(abs(diffs - median(diffs)));
    sigma_est = mad_val / 0.6745;
    threshold = 3 * sigma_est;

    Z_corr = Z;
    cum_offset = 0;
    for k = 2:nrows
        if abs(diffs(k-1)) > threshold
            cum_offset = cum_offset - diffs(k-1);
        end
        Z_corr(k, :) = Z(k, :) + cum_offset;
    end
end

% Sliding-window std dev using var = <x^2> - <x>^2
function V = local_std_map(Z, win)
    k    = ones(win, win) / win^2;
    mu   = conv2(Z,    k, 'same');
    mu2  = conv2(Z.^2, k, 'same');
    V    = sqrt(max(mu2 - mu.^2, 0));
end

% colormap
function cmap = bwr_cmap(n)
    half = floor(n / 2);
    blue_to_white = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    white_to_red  = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
    cmap = [blue_to_white; white_to_red];
end

% color hex conversion
function rgb = hex2rgb(hex)
    hex = strrep(hex, '#', '');
    rgb = double([hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]) / 255;
end
