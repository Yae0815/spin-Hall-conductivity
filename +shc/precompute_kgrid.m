function cache = precompute_kgrid(params)
% shc.precompute_kgrid
% Precompute E, velocity/spin matrix elements, and AB/BA denominators.

    arguments
        params struct
    end

    Nk   = params.Nk;
    eta  = params.eta;
    hbar = params.hbar;
    alpha = lower(params.alpha);
    beta  = lower(params.beta);
    gamma = lower(params.gamma);

    build = shc.make_builders(params);
    Norb  = build.Norb;
    assert(Norb==4, 'This precompute is optimized for Norb=4.');

    % k-grid (centered)
    base = (0:Nk-1)/Nk - 0.5;
    [KX,KY,KZ] = ndgrid(2*pi*base, 2*pi*base, 2*pi*base);
    K = numel(KX);

    E   = zeros(Norb, K);
    Xab = zeros(Norb, Norb, K);
    Xba = zeros(Norb, Norb, K);
    Den = zeros(Norb, Norb, K);
    w   = ones(K,1) * ( (2*pi)^3 / Nk^3 );  % simple equal weight

    Sg = shc.util_pick_spin_op(build, gamma);

    parfor t = 1:K
        k = [KX(t), KY(t), KZ(t)];
        Hk  = build.H(k);
        [U,D] = eig((Hk+Hk')/2);            % hermitize (numerical safe)
        evals = real(diag(D));
        E(:,t) = evals;

        % velocities
        Vk_a = build.vel(alpha, k, hbar);
        Vk_b = build.vel(beta , k, hbar);

        % matrix elements in band basis
        Va = U' * Vk_a * U;
        Vb = U' * Vk_b * U;
        S  = U' * Sg  * U;

        % AB-BA (spin-velocity commutators in Kubo)
        Xab(:,:,t) = (S*Va)*Vb;
        Xba(:,:,t) = (S*Vb)*Va;

        % denominators (E_n - E_m)^2 + eta^2  for n,m
        En = evals(:);
        dE = En - En.';
        Den(:,:,t) = dE.^2 + (eta^2)*eye(Norb);
    end

    cache = struct('E',E, 'Xab',Xab, 'Xba',Xba, 'Den',Den, ...
                   'w',w, 'Nk',Nk, 'eta',eta, 'alpha',alpha, 'beta',beta, ...
                   'gamma',gamma, 'build',build);
end
