%% ===== sweep_Nk_shifted.m =====
clear;
load ../ftn58sparse.mat

% 固定參數
base.ftn58 = ftn58sparse;
base.Ef    = 0.0;
base.eta   = 1e-4;
base.electronic_charge = 1.0;
base.hbar  = 1.0;


Nk_list = [40 60 80 100 130 150 180];

% 建議的 5 個位移（單位：格點），效果好且成本低
shifts = [
  0 0 0;  0.5 0 0;  0 0.5 0;  0 0 0.5;
  0.5 0.5 0;  0.5 0 0.5;  0 0.5 0.5;  0.5 0.5 0.5
];

% 可選：啟動平行（spin_hall_main 內部已有 parfor；外層不用再 parfor）
try
    if isempty(gcp('nocreate')), parpool('threads'); end
end

% 儲存容器
sx_xz_mean = zeros(size(Nk_list));
sx_xz_sem  = zeros(size(Nk_list));   % standard error of mean
sy_yz_mean = zeros(size(Nk_list));
sy_yz_sem  = zeros(size(Nk_list));
timesec    = zeros(size(Nk_list));

fprintf('Nk  |  <σ^{s_x}_{xz}> ± SEM         <σ^{s_y}_{yz}> ± SEM         time(s)\n');
fprintf('----+---------------------------------------------------------------\n');

for t = 1:numel(Nk_list)
    params = base;
    params.Nk = Nk_list(t);

    vals_x = zeros(size(shifts,1),1);
    vals_y = zeros(size(shifts,1),1);

    t0 = tic;
    for si = 1:size(shifts,1)
        params.shift = shifts(si,:);

        % σ^{s_x}_{xz}
        params.alpha = 'x'; params.beta = 'z'; params.gamma = 'x';
        res = spin_hall_main(params);
        vals_x(si) = res.sigma_sab_gamma;

        % σ^{s_y}_{yz}
        params.alpha = 'y'; params.beta = 'z'; params.gamma = 'y';
        res = spin_hall_main(params);
        vals_y(si) = res.sigma_sab_gamma;
    end
    timesec(t) = toc(t0);

    % 平均與標準誤
    sx_xz_mean(t) = mean(vals_x);
    sy_yz_mean(t) = mean(vals_y);
    sx_xz_sem(t)  = std(vals_x,0,1) / sqrt(numel(vals_x));
    sy_yz_sem(t)  = std(vals_y,0,1) / sqrt(numel(vals_y));

    fprintf('%3d |  % .6e ± %.1e    % .6e ± %.1e    %7.2f\n', ...
        params.Nk, sx_xz_mean(t), sx_xz_sem(t), sy_yz_mean(t), sy_yz_sem(t), timesec(t));
end

save('spin_hall_convergence_shift_avg.mat', ...
     'Nk_list','sx_xz_mean','sx_xz_sem','sy_yz_mean','sy_yz_sem','timesec','shifts','base');