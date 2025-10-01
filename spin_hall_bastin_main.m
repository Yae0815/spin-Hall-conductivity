function out = spin_hall_bastin_main(params)
% Bastin (f_n - f_m) 主流程
    if ~isfield(params,'method'), params.method = 'bastin'; end
    params.method = 'bastin';

    if ~isfield(params,'Ef'),      params.Ef      = 0.0; end
    if ~isfield(params,'T'),       params.T       = 0;   end
    if ~isfield(params,'mu_list'), params.mu_list = 0;   end

    cache  = shc.precompute_kgrid(params);
    mu_vec = params.mu_list(:).';
    sigma  = zeros(size(mu_vec));
    for i = 1:numel(mu_vec)
        sigma(i) = shc.eval_sigma(cache, mu_vec(i), params.Ef, params.T, 'bastin');
    end

    H0   = cache.build.H([0,0,0]);
    Sg   = shc.util_pick_spin_op(cache.build, params.gamma);
    comm = Sg*H0 - H0*Sg;
    comm_norm = norm(comm, 'fro');

    out = struct('mu',mu_vec,'sigma',sigma,'params',params,'cache',cache, ...
                 'sigma_sab_gamma', sigma(1), 'commutator_norm_Gamma', comm_norm);
end
