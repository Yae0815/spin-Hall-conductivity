% plot_all_results.m
clear; clf;
files = dir('result_mu_scan_*.mat');  % 找出所有結果檔
hold on;

for k = 1:numel(files)
    data = load(files(k).name, 'result');
    r = data.result;
    
    plot(r.mu_grid, r.sig, 'LineWidth', 1.6, ...
         'DisplayName', r.legend);   % 用你手動輸入的名稱
end

% 統一標籤
xlabel('$\mu\ \mathrm{(eV)}$', 'Interpreter','latex','FontSize',23);
ylabel('$\sigma^{s_x}_{xz}\ \mathrm{(lattice\ units)}$', ...
       'Interpreter','latex','FontSize',23);
title('$\mathrm{Spin\ Hall\ vs.\ Chemical\ Potential}$', ...
      'Interpreter','latex','FontSize',23);

ax = gca; ax.FontSize = 18;
grid on; box on;
legend('Interpreter','latex','FontSize',15, 'Location','best');

