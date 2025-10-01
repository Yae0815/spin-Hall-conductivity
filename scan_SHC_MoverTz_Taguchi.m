%% scan_SHC_MoverTz_Taguchi.m
% 掃描 M/tz 對 σ^{s_z}_{xy} 的影響，並對多組 txy/tz 繪製曲線（Taguchi 2020 DSM）
% 需求：
%   - 已有：spin_hall_main.m（你的現成內核）
%   - 已有：make_builders_taguchi.m（上回給你的 Taguchi builder）
% 使用：
%   直接 Run 本檔；輸出 .mat 與 .png

clear; clc;

%% ===== User controls =====
% --- Lattice model ratios（建議覆現區間，含 TDSM/DSM 轉換） ---
txy_over_tz_list = [0.5, 1.0, 1.5];    % 多條曲線
M_over_tz_grid   = linspace(0.0, 2.0, 101);

% --- 固定的比例（對照 Taguchi 常用設定） ---
eta_unit      = 0.89;          % 整體能量尺度（eV）
tz_over_eta   = -3.4;          % tz/eta
beta_over_tz  = 0.67;          % β/tz
gamma_over_tz = 0.335;         % γ/tz (= (1/2)*β/tz)

% --- 數值積分與展寬 ---
params.Nk   = 101;              % odd; 增大可更平滑
params.eta  = 1e-4;            % Kubo broadening (eV)
params.Ef   = 0.0;
params.electronic_charge = 1.0;
params.hbar = 1.0;

% --- 只算這一個張量元素 ---
params.builder = 'taguchi2020';
params.alpha = 'x';
params.beta  = 'y';
params.gamma = 'z';

%% ===== Derived constants =====
tz = tz_over_eta * eta_unit;   % eV
beta  = beta_over_tz  * tz;
gamma = gamma_over_tz * tz;

nM   = numel(M_over_tz_grid);
nRat = numel(txy_over_tz_list);

SHC  = zeros(nM, nRat);        % 每欄對應一個 txy/tz；列是 M/tz
comm_at_Gamma = zeros(nRat,1); % 只記第一個 M 作 sanity check

fprintf('Scan σ^{s_z}_{xy} vs M/tz for several txy/tz ... Nk=%d, eta=%.1e\n', params.Nk, params.eta);

t_start = tic;
for ir = 1:nRat
    txy_over_tz = txy_over_tz_list(ir);
    txy = txy_over_tz * tz;

    % 先做一次（用中間的 M）建立 baseline，並記錄 Γ 的 [Sz,H] 大小
    midM = M_over_tz_grid( ceil(nM/2) ) * tz;
    params.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',midM,'beta',beta,'gamma',gamma);
    res0 = spin_hall_main(params);
    comm_at_Gamma(ir) = res0.commutator_norm_Gamma;

    % 掃描 M/tz
    par_SHC = zeros(nM,1);
    parfor iM = 1:nM
        Mi = M_over_tz_grid(iM) * tz;
        p_local = params;
        p_local.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',Mi,'beta',beta,'gamma',gamma);

        r = spin_hall_main(p_local);          % 只算 σ^{s_z}_{xy}
        par_SHC(iM) = r.sigma_sab_gamma;      % lattice units
    end
    SHC(:,ir) = par_SHC;

    fprintf(' done: txy/tz = %.3f  (||[S_z,H(Γ)]||_F ~ %.2e)\n', txy_over_tz, comm_at_Gamma(ir));
end
elapsed = toc(t_start);

%% ===== Plot =====
figure('Color','w');
co = [1 0 0; 0 0 1; 0.5 0 0.5];   % 紅、藍、紫 (RGB)
set(gca,'ColorOrder',co,'NextPlot','replacechildren');
hold on;
for ir = 1:nRat
    plot(M_over_tz_grid, SHC(:,ir), 'LineWidth', 1.8);
end
grid on;
xlabel('M/t_z','FontSize',12);
ylabel('\sigma^{s_z}_{xy} (lattice units)','FontSize',12);
legstr = arrayfun(@(x) sprintf('t_{xy}/t_z = %.3f', x), txy_over_tz_list, 'UniformOutput', false);
legend(legstr, 'Location','best', 'Interpreter','tex');
title(sprintf('Spin Hall: \\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV)', params.Nk, params.eta));

% 存檔
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
png_name = sprintf('SHC_MoverTz_scan_%s.png', ts);
mat_name = sprintf('SHC_MoverTz_scan_%s.mat', ts);

exportgraphics(gcf, png_name, 'Resolution', 200);
meta = struct( ...
    'Nk', params.Nk, 'broadening_eV', params.eta, 'Ef', params.Ef, ...
    'eta_unit', eta_unit, 'tz_over_eta', tz_over_eta, ...
    'beta_over_tz', beta_over_tz, 'gamma_over_tz', gamma_over_tz, ...
    'txy_over_tz_list', txy_over_tz_list, 'M_over_tz_grid', M_over_tz_grid, ...
    'elapsed_sec', elapsed, 'note', 'sigma^{s_z}_{xy} only; Taguchi lattice model' ...
);
save(mat_name, 'SHC', 'meta', 'comm_at_Gamma');

fprintf('\nSaved: %s  and  %s\n', mat_name, png_name);

%% ===== Plot (Y data inverted: -σ) =====
figure('Color','w');
co = [1 0 0; 0 0 1; 0.5 0 0.5];   % 紅、藍、紫 (RGB)
set(gca,'ColorOrder',co,'NextPlot','replacechildren');
hold on;
for ir = 1:nRat
    plot(M_over_tz_grid, -SHC(:,ir), 'LineWidth', 1.8);  % 注意這裡是負號
end
grid on;
xlabel('M/t_z','FontSize',12);
ylabel('-\sigma^{s_z}_{xy} (lattice units)','FontSize',12);
ylim([-0.07 0.01])
legstr = arrayfun(@(x) sprintf('t_{xy}/t_z = %.3f', x), txy_over_tz_list, 'UniformOutput', false);
legend(legstr, 'Location','best', 'Interpreter','tex');
title(sprintf('Spin Hall: -\\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV)', params.Nk, params.eta));

png_name_inv = sprintf('SHC_MoverTz_scan_yINV_%s.png', ts);
exportgraphics(gcf, png_name_inv, 'Resolution', 200);
fprintf('Saved inverted plot: %s\n', png_name_inv);
