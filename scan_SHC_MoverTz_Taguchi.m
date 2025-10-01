%% scan_SHC_MoverTz_Taguchi.m
% 掃描 M/tz 對 σ^{s_z}_{xy} 的影響，並對多組 txy/tz 繪製曲線（Taguchi 2020 DSM）
% - 入口使用 spin_hall_main(params)
% - 單位：可選 'lattice' 或 'SI'（對齊論文圖3(b)：(ħ/e)(Ω·m)^-1）
% - builder: 'taguchi2020'（需 make_builders_taguchi.m）

clear; clc;

%% ===== User controls =====
% --- Lattice model ratios（建議覆現區間，含 TDSM/DSM 轉換） ---
txy_over_tz_list = [0.5, 1.0, 1.5];      % 多條曲線
M_over_tz_grid   = linspace(0.0, 2.0, 21);

% --- 固定的比例（對照 Taguchi 常用設定） ---
eta_unit      = 0.89;          % 整體能量尺度（eV）
tz_over_eta   = -3.4;          % tz/eta
beta_over_tz  = 0.67;          % β/tz
gamma_over_tz = 0.335;         % γ/tz (= (1/2)*β/tz)

% --- 數值積分與展寬 ---
params.Nk   = 41;              % odd; 增大可更平滑
params.eta  = 1e-4;            % Kubo broadening (eV)
params.Ef   = 0.0;
params.electronic_charge = 1.0;
params.hbar = 1.0;

% --- 只算這一個張量元素 ---
params.builder = 'taguchi2020';
params.alpha = 'x';
params.beta  = 'y';
params.gamma = 'z';

% --- Kubo 實作分支（若 spin_hall_main 有接 method 參數會使用） ---
params.method = 'weighted';     % 'weighted' 或 'bastin'

% --- 輸出單位（'lattice' 或 'SI'）---
out_units = 'SI';               % 對齊論文圖3(b)：(ħ/e)(Ω·m)^-1
a_angstrom = 1.0;               % 晶格常數 a（Å）；請依材料調整

%% ===== Derived constants =====
tz = tz_over_eta * eta_unit;   % eV
beta  = beta_over_tz  * tz;
gamma = gamma_over_tz * tz;

nM   = numel(M_over_tz_grid);
nRat = numel(txy_over_tz_list);

SHC_lattice = zeros(nM, nRat);      % lattice units（無因次 S）
comm_at_Gamma = zeros(nRat,1);      % 只記第一個 M 作 sanity check

% 轉換係數：S -> 數值（圖3(b)的 y 軸數字），軸標為 [(ħ/e)(Ω·m)^-1]
% scale = (e^2/h)/a
e2_over_h_S = 3.874045e-5;          % Siemens
a_meter = a_angstrom * 1e-10;       % m
scale_SI = e2_over_h_S / a_meter;   % (Ω·m)^-1

fprintf('Scan σ^{s_z}_{xy} vs M/tz for several txy/tz ... Nk=%d, eta=%.1e, method=%s\n', ...
    params.Nk, params.eta, params.method);

t_start = tic;
for ir = 1:nRat
    txy_over_tz = txy_over_tz_list(ir);
    txy = txy_over_tz * tz;

    % 先做一次（用中間的 M）建立 baseline，並記錄 Γ 的 [Sz,H] 大小
    midM = M_over_tz_grid( ceil(nM/2) ) * tz;
    params.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',midM,'beta',beta,'gamma',gamma);
    res0 = spin_hall_main(params);
    if isfield(res0,'commutator_norm_Gamma')
        comm_at_Gamma(ir) = res0.commutator_norm_Gamma;
    else
        comm_at_Gamma(ir) = NaN;
    end

    % 掃描 M/tz
    par_SHC = zeros(nM,1);
    parfor iM = 1:nM
        Mi = M_over_tz_grid(iM) * tz;
        p_local = params;
        p_local.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',Mi,'beta',beta,'gamma',gamma);

        r = spin_hall_main(p_local);              % 只算 σ^{s_z}_{xy}
        % r.sigma_sab_gamma 應為 lattice units 的無因次 S
        par_SHC(iM) = r.sigma_sab_gamma;
    end
    SHC_lattice(:,ir) = par_SHC;

    fprintf(' done: txy/tz = %.3f  (||[S_z,H(Γ)]||_F ~ %.2e)\n', txy_over_tz, comm_at_Gamma(ir));
end
elapsed = toc(t_start);

%% ===== Unit conversion for plotting =====
switch lower(out_units)
    case 'si'
        Y = SHC_lattice .* scale_SI;  % 數值直接是圖3(b)的 y 值；軸標標註 (ħ/e)(Ω·m)^-1
        ylab = '\sigma^{s_z}_{xy}  [ (ħ/e)\,(\Omega\cdot m)^{-1} ]';
        units_tag = 'SI';
    case 'lattice'
        Y = SHC_lattice;
        ylab = '\sigma^{s_z}_{xy} (lattice units)';
        units_tag = 'lattice';
    otherwise
        error('Unknown out_units: %s', out_units);
end

%% ===== Plot =====
figure('Color','w');
hold on;
for ir = 1:nRat
    plot(M_over_tz_grid, Y(:,ir), 'LineWidth', 1.8);
end
grid on;
xlabel('M/t_z','FontSize',12);
ylabel(ylab,'FontSize',12);
legstr = arrayfun(@(x) sprintf('t_{xy}/t_z = %.3f', x), txy_over_tz_list, 'UniformOutput', false);
legend(legstr, 'Location','best', 'Interpreter','tex');
title(sprintf('Spin Hall: \\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV, %s)', ...
      params.Nk, params.eta, params.method));

%% ===== Save =====
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
png_name = sprintf('SHC_MoverTz_scan_%s_%s.png', ts, units_tag);
mat_name = sprintf('SHC_MoverTz_scan_%s_%s.mat', ts, units_tag);

exportgraphics(gcf, png_name, 'Resolution', 220);

meta = struct( ...
    'Nk', params.Nk, 'broadening_eV', params.eta, 'Ef', params.Ef, ...
    'eta_unit', eta_unit, 'tz_over_eta', tz_over_eta, ...
    'beta_over_tz', beta_over_tz, 'gamma_over_tz', gamma_over_tz, ...
    'txy_over_tz_list', txy_over_tz_list, 'M_over_tz_grid', M_over_tz_grid, ...
    'elapsed_sec', elapsed, 'note', 'sigma^{s_z}_{xy} only; Taguchi lattice model', ...
    'units', out_units, 'a_angstrom', a_angstrom, ...
    'scale_SI_per_lattice', scale_SI, ...
    'comm_at_Gamma', comm_at_Gamma, ...
    'method', params.method ...
);
save(mat_name, 'SHC_lattice', 'Y', 'meta');

fprintf('\nSaved: %s  and  %s\n', mat_name, png_name);

%% ===== Plot (-Y) =====
Y_neg = -Y;

figure('Color','w');
hold on;
for ir = 1:nRat
    plot(M_over_tz_grid, Y_neg(:,ir), 'LineWidth', 1.8);
end
grid on;
xlabel('M/t_z','FontSize',12);
% y 標籤沿用同一單位，但在 title 註記「-Y」
ylabel(ylab,'FontSize',12);
legend(legstr, 'Location','best', 'Interpreter','tex');
title(sprintf('Spin Hall: -\\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV, %s)', ...
      params.Nk, params.eta, params.method));

%% ===== Save (-Y) =====
png_name_neg = sprintf('SHC_MoverTz_scan_%s_%s_neg.png', ts, units_tag);
exportgraphics(gcf, png_name_neg, 'Resolution', 220);

% （可選）把 Y_neg 一起存進 .mat 方便之後重畫
if exist('mat_name','var') && isfile(mat_name)
    S = whos('-file', mat_name);
    tosave = ~ismember('Y_neg', {S.name});
    if tosave
        save(mat_name, 'Y_neg', '-append');
    end
else
    % 若你想單獨另存一份含 Y_neg 的 .mat（不覆蓋原檔），用下面這行：
    % save(strrep(mat_name, '.mat', '_withNeg.mat'), 'SHC_lattice','Y','Y_neg','meta');
end

fprintf('Saved (negative plot): %s\n', png_name_neg);