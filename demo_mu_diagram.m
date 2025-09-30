% demo_mu_diagram.m  (precompute + μ-scan; DSM sanity check)
clear; close all;
tic;  % start clock

% ===== add path =====
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, 'shc_kubo')));  % add +shc

% ===== params: Spin Hall tensor indices (strings) =====
params.alpha              = 'x';   % v_α
params.beta               = 'y';   % v_β
params.gamma              = 'z';   % S_γ

% ===== grid / numerics =====
params.Nk                 = 41;
params.eta                = 5e-3;          % eV (Kubo broadening)
params.electronic_charge  = 1;             % lattice units
params.hbar               = 1;             % lattice units
params.shift              = [0 0 0];       % e.g. [0.5 0.5 0.5]/Nk

Ef                        = 3.0;           % reference Fermi energy

% ===== model (Taguchi-like DSM) — 注意：模型用 beta_cpl / gamma_cpl =====
params.eta_band           = 0.1;
params.txy                = 1.5;
params.tz                 = 1.0;
params.beta_cpl           = 0.0;           % <-- 原先叫 beta 的模型耦合，為避免撞名改為 beta_cpl
params.gamma_cpl          = 0.0;           % <-- 原先叫 gamma 的模型耦合，改為 gamma_cpl
params.M                  = 1.0;

% ===== 1) precompute once =====
cache = shc.precompute_kgrid(params);


% ===== 2) evaluate for many μ (and T if wanted) =====
mu_grid = linspace(-4, 4, 41);
T       = 10;                  % K (e.g., set 300 for finite-T test)
method  = 'weighted';          % or 'bastin'

sig = zeros(numel(mu_grid),1);
for i = 1:numel(mu_grid)
    out    = shc.eval_sigma(cache, mu_grid(i), Ef, T, method);
    sig(i) = out.sigma;
end
toc;  % stop clock

% ===== 3) plot =====
figure('Color','w'); plot(mu_grid, sig, 'LineWidth',1.6);
grid on; box on;

lbl_interpreter = 'latex';
unit_note       = '\mathrm{(lattice\ units)}';

a = char(lower(string(params.alpha)));
b = char(lower(string(params.beta)));
g = char(lower(string(params.gamma)));

yl_core = sprintf('\\sigma^{s_{%s}}_{%s%s}', g, a, b);
yl = ['$' yl_core '\; ' unit_note '$'];
xl = '$\mu\ \mathrm{(eV)}$';

xlabel(xl, 'Interpreter', lbl_interpreter, 'FontSize', 23);
ylabel(yl, 'Interpreter', lbl_interpreter, 'FontSize', 23);
title(sprintf('$\\mathrm{Spin\\ Hall\\ vs.\\ Chemical\\ Potential}\\ (T=%.1f\\,\\mathrm{K})$', T), ...
      'Interpreter','latex', 'FontSize', 23);

ax = gca; ax.FontSize = 18;

fprintf('Precompute done: Nk=%d, eta=%.3g eV. Eval method=%s, T=%.1f K\n', ...
    params.Nk, params.eta, method, T);
beep;

% ===== 4) save results =====
result.mu_grid = mu_grid;
result.sig     = sig;
result.params  = params;
result.T       = T;
result.method  = method;
result.Ef      = Ef;

result.legend  = input('Enter legend label for this run: ','s');

timestamp = string(datetime("now","Format","yyyyMMdd_HHmmss"));
fname = sprintf('result_mu_scan_%s.mat', timestamp);
save(fname, 'result');
fprintf('Saved results to %s\n', fname);
