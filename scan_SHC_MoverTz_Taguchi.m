%% scan_SHC_MoverTz_Taguchi.m
clear; clc;

%% ===== Path setup =====
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));

%% ===== User controls =====
method_choice    = 'bastin';           % 'weighted' or 'bastin'
txy_over_tz_list = [0.5, 1.0, 1.5];
M_over_tz_grid   = linspace(0.0, 2.0, 101);

eta_unit      = 0.89;          % eV
tz_over_eta   = -3.4;
beta_over_tz  = 0.67;
gamma_over_tz = 0.335;

params.Nk   = 51;
params.eta  = 1e-4;
params.Ef   = 0.0;
params.electronic_charge = 1.0;
params.hbar = 1.0;

params.builder = 'taguchi2020';
params.alpha   = 'x';
params.beta    = 'y';
params.gamma   = 'z';
params.method  = method_choice;

%% ===== Derived constants =====
tz = tz_over_eta * eta_unit;
beta  = beta_over_tz  * tz;
gamma = gamma_over_tz * tz;

nM   = numel(M_over_tz_grid);
nRat = numel(txy_over_tz_list);

SHC  = zeros(nM, nRat);
comm_at_Gamma = zeros(nRat,1);

switch lower(method_choice)
    case 'weighted', runner = @spin_hall_main;
    case 'bastin',   runner = @spin_hall_bastin_main;
    otherwise, error('Unknown method: %s', method_choice);
end

fprintf('Scan σ^{s_z}_{xy} vs M/tz ... method=%s, Nk=%d, eta=%.1e\n', ...
    method_choice, params.Nk, params.eta);

t_start = tic;
for ir = 1:nRat
    txy_over_tz = txy_over_tz_list(ir);
    txy = txy_over_tz * tz;

    midM = M_over_tz_grid( ceil(nM/2) ) * tz;
    params.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',midM,'beta',beta,'gamma',gamma);
    res0 = runner(params);
    comm_at_Gamma(ir) = res0.commutator_norm_Gamma;

    par_SHC = zeros(nM,1);
    parfor iM = 1:nM
        Mi = M_over_tz_grid(iM) * tz;
        p_local = params;
        p_local.model = struct('eta',eta_unit,'txy',txy,'tz',tz,'M',Mi,'beta',beta,'gamma',gamma);

        r = runner(p_local);
        par_SHC(iM) = r.sigma_sab_gamma;
    end
    SHC(:,ir) = par_SHC;

    fprintf(' done: txy/tz = %.3f  (||[S_z,H(Γ)]||_F ~ %.2e)\n', txy_over_tz, comm_at_Gamma(ir));
end
elapsed = toc(t_start);

%% ===== Plot =====
figure('Color','w');
co = [1 0 0; 0 0 1; 0.5 0 0.5];
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
title(sprintf('Spin Hall (%s): \\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV)', ...
      upper_first(method_choice), params.Nk, params.eta));

ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
png_name = sprintf('SHC_MoverTz_scan_%s_%s.png', method_choice, ts);
mat_name = sprintf('SHC_MoverTz_scan_%s_%s.mat', method_choice, ts);

exportgraphics(gcf, png_name, 'Resolution', 200);
meta = struct('method',method_choice,'Nk',params.Nk,'broadening_eV',params.eta,'Ef',params.Ef, ...
              'eta_unit',eta_unit,'tz_over_eta',tz_over_eta,'beta_over_tz',beta_over_tz,'gamma_over_tz',gamma_over_tz, ...
              'txy_over_tz_list',txy_over_tz_list,'M_over_tz_grid',M_over_tz_grid,'elapsed_sec',elapsed, ...
              'note','sigma^{s_z}_{xy} only; Taguchi lattice model');
save(mat_name, 'SHC', 'meta', 'comm_at_Gamma');
fprintf('\nSaved: %s  and  %s\n', mat_name, png_name);

%% ===== Plot (Y inverted) =====
figure('Color','w');
set(gca,'ColorOrder',co,'NextPlot','replacechildren');
hold on;
for ir = 1:nRat
    plot(M_over_tz_grid, -SHC(:,ir), 'LineWidth', 1.8);
end
grid on;
xlabel('M/t_z','FontSize',12);
ylabel('-\sigma^{s_z}_{xy} (lattice units)','FontSize',12);
ylim([-0.07 0.01])
legend(legstr, 'Location','best', 'Interpreter','tex');
title(sprintf('Spin Hall (%s): -\\sigma^{s_z}_{xy} vs M/t_z  (Nk=%d, \\eta=%.0e eV)', ...
      upper_first(method_choice), params.Nk, params.eta));

png_name_inv = sprintf('SHC_MoverTz_scan_yINV_%s_%s.png', method_choice, ts);
exportgraphics(gcf, png_name_inv, 'Resolution', 200);
fprintf('Saved inverted plot: %s\n', png_name_inv);

function s = upper_first(str)
    if isempty(str), s = str; return; end
    s = lower(string(str));
    s = extractBetween(s,1,1).upper + extractAfter(s,1);
    s = char(s);
end
