function scan = scan_mu(params, mu_grid)
% shc.scan_mu
% Sweep chemical potential μ and record σ^{s_γ}_{αβ}(μ).
% NOTE: spin_hall_main already uses parfor over kx,
% so we keep a plain for-loop here to avoid nested parallelism.

    sig = zeros(numel(mu_grid),1);

    for i = 1:numel(mu_grid)
        p = params;
        p.Ef = mu_grid(i);
        out = shc.spin_hall_main(p);
        sig(i) = out.sigma_sab_gamma;
    end

    scan.mu    = mu_grid(:);
    scan.sigma = sig(:);
end
