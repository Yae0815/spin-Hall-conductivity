function out = spin_hall_bastin_main(params)
% spin_hall_bastin_main
% Run Spin Hall conductivity using the Bastin formula (f_n - f_m).
%
% Usage:
%   params = struct;
%   % ... fill your params (Nk, eta, alpha,beta,gamma, Ef, mu_list, T, etc.)
%   out = spin_hall_bastin_main(params);

    arguments
        params struct
    end

    % ---- defaults (only set if missing) ----
    if ~isfield(params, 'method'), params.method = 'bastin'; end
    params.method = 'bastin';   % force Bastin

    if ~isfield(params, 'Ef'),      params.Ef      = 0.0;   end
    if ~isfield(params, 'T'),       params.T       = 0;     end
    if ~isfield(params, 'mu_list'), params.mu_list = 0;     end

    % ===== precompute on k-grid =====
    cache = shc.precompute_kgrid(params);   % expects your existing +shc code

    % ===== sweep μ =====
    mu_list = params.mu_list(:).';
    sigma   = zeros(size(mu_list));
    for i = 1:numel(mu_list)
        sigma(i) = shc.eval_sigma(cache, mu_list(i), params.Ef, params.T, 'bastin');
    end

    % ===== pack & (optionally) plot =====
    out = struct('mu', mu_list, 'sigma', sigma, 'params', params, 'cache', cache);

    if ~isfield(params,'no_plot') || ~params.no_plot
        figure; plot(mu_list, sigma, 'LineWidth', 2);
        grid on; xlabel('\mu (eV)'); ylabel('\sigma^{s_\gamma}_{\alpha\beta}');
        title('Spin Hall (Bastin)');
    end
end
