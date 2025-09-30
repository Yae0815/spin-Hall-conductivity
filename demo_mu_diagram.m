% demo_mu_diagram.m  (precompute + μ-scan)
clear;
tic; % start clock
ticBytes(gcp);
% ===== add path =====
%this_dir = fileparts(mfilename('fullpath'));
%addpath(genpath(fullfile(this_dir, '..'))); % add shc_kubo/+shc

% ===== load model =====
%Sanity Check
%ftn = load('../ftn58sparse.mat');
%if isfield(ftn,'ftn58sparse'), ftn = ftn.ftn58sparse; end

% ===== params =====
%Sanity Check
%params.ftn58              = ftn;
params.Nk                 = 21;
params.eta                = 5e-3;           % eV (fixed in cache; change => re-precompute)
params.electronic_charge  = 1;              % lattice units
params.hbar               = 1;              % lattice units
params.alpha              = 'x';
params.beta               = 'y';
params.gamma              = 'z';
params.shift              = [0 0 0];        % can test [0.5 0.5 0.5]/Nk

Ef                        = 3.0;

% ===== 0) DSM =====
params.eta_band = 0.1;   % <-- 用於 a1,a2 的 η（與 Kubo 展寬不同名）
params.txy = 1.5;
params.tz = 1;
params.beta = 0;
params.gamma = 0;
params.M = 1;
% ===== 1) precompute once =====
cache = shc.precompute_kgrid(params);

% ===== 2) evaluate for many μ (and T if wanted) =====
mu_grid = linspace(-4, 4, 101);
T       = 10;                  % change to e.g. 300 (K) to test finite-T
method  = 'weighted';         % or 'bastin'

sig = zeros(numel(mu_grid),1);
for i = 1:numel(mu_grid)
    out = shc.eval_sigma(cache, mu_grid(i),Ef, T, method);
    sig(i) = out.sigma;
end

toc; % stop clock
tocBytes(gcp)
% ===== 3) plot =====
figure('Color','w'); plot(mu_grid, sig, 'LineWidth',1.6);

% labels (純 lattice units，去掉 e/ħ)
lbl_interpreter = 'latex';                      % 'latex' 或 'tex'
unit_note       = '\mathrm{(lattice\ units)}';

a = lower(params.alpha); b = lower(params.beta); g = lower(params.gamma);
switch lower(lbl_interpreter)
    case 'latex'
        yl_core = sprintf('\\sigma^{s_{%s}}_{%s%s}', g, a, b);
        yl = ['$' yl_core '\; ' unit_note '$'];
        xl = '$\mu\ \mathrm{(eV)}$';
    otherwise
        yl = sprintf('\\sigma^{s_{%s}}_{%s%s} %s', g, a, b, unit_note);
        xl = '\mu (eV)';
end


xlabel(xl, 'Interpreter', lbl_interpreter, 'FontSize', 23);
ylabel(yl, 'Interpreter', lbl_interpreter, 'FontSize', 23);
title(sprintf('$\\mathrm{Spin\\ Hall\\ vs.\\ Chemical\\ Potential}\\ (T=%.1f\\,\\mathrm{K})$', T), ...
      'Interpreter','latex', 'FontSize', 23);

ax = gca;
ax.FontSize = 18;   
grid on; box on;

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

% 手動輸入 legend 名稱
result.legend  = input('Enter legend label for this run: ','s');

% 自動產生檔名
timestamp = string(datetime("now","Format","yyyyMMdd_HHmmss"));
fname = sprintf('result_mu_scan_%s.mat', timestamp);
save(fname, 'result');
fprintf('Saved results to %s\n', fname);
