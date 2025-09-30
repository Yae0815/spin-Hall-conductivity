function cache = precompute_kgrid(params)
% shc.precompute_kgrid
% One-time preprocessing for Spin Hall (Kubo AB-BA) with N=4.
% Stores, for each k, E, Xab, Xba, and denom matrix.
% Changing eta/alpha/beta/gamma/Nk/shift requires re-precompute.

    % ---- unpack indices (robust to string/char) ----
    alpha = char(lower(string(params.alpha)));
    beta  = char(lower(string(params.beta)));
    gamma = char(lower(string(params.gamma)));
    ok = @(s) ischar(s) && ismember(s, {'x','y','z'});
    assert(ok(alpha) && ok(beta) && ok(gamma), ...
        'params.alpha/beta/gamma must be ''x'',''y'',''z'' strings.');

    % ---- numerics ----
    Nk    = params.Nk;
    eta   = params.eta;
    hbar  = params.hbar;

    % ---- builders (DSM sanity model) ----
    build = shc.shc_builders_dsm_taguchi(params);
    Norb  = build.Norb;
    assert(Norb==4, 'This precompute is optimized for Norb=4.');

    % ---- spin operator from builder ----
    if isfield(build,'S') && isfield(build.S, gamma)
        Sg = build.S.(gamma);
    elseif all(isfield(build, {'Sx','Sy','Sz'}))
        switch gamma
            case 'x', Sg = build.Sx;
            case 'y', Sg = build.Sy;
            case 'z', Sg = build.Sz;
        end
    else
        error('Spin operators not found in builder. Expect build.S.(x|y|z) or build.Sx/Sy/Sz.');
    end

    % ---- k grid ----
    base  = (0:Nk-1)/Nk - 0.5;
    shift = [0 0 0];
    if isfield(params,'shift') && ~isempty(params.shift), shift = params.shift; end
    kx_list = base + shift(1)/Nk;
    ky_list = base + shift(2)/Nk;
    kz_list = base + shift(3)/Nk;
    w = 1.0 / (Nk^3);

    % ---- handles ----
    Hk  = build.H;
    dHx = build.dHdkx; dHy = build.dHdky; dHz = build.dHdkz;

    % ---- allocate (4x4xK) & (4xK) ----
    K = Nk*Nk*Nk;
    Xab = zeros(4,4,K);   % 2*imag(J_A .* v_B.')
    Xba = zeros(4,4,K);
    Den = zeros(4,4,K);   % (E_m - E_n)^2 + eta^2, diag = Inf
    Eall= zeros(4,K);     % eigenvalues

    % ---- flatten k loop for parfor ----
    [IX,IY,IZ] = ndgrid(1:Nk, 1:Nk, 1:Nk);
    ix = IX(:); iy = IY(:); iz = IZ(:);

    parfor t = 1:K
        kx = kx_list(ix(t)); ky = ky_list(iy(t)); kz = kz_list(iz(t));

        % H and velocities (via unified picker)
        H   = Hk(kx,ky,kz);
        vA  = shc.pick_velocity(dHx,dHy,dHz, alpha, kx,ky,kz, hbar);
        vB  = shc.pick_velocity(dHx,dHy,dHz, beta,  kx,ky,kz, hbar);

        % hermitize just in case
        vA  = (vA + vA')/2; vB = (vB + vB')/2;

        % spin current operators
        JAs = 0.5*(Sg*vA + vA*Sg);
        JBs = 0.5*(Sg*vB + vB*Sg);
        JAs = (JAs + JAs')/2; JBs = (JBs + JBs')/2;

        % diagonalize H
        [U,D] = eig(full(H));
        E = real(diag(D)); [E,idx] = sort(E,'ascend'); U = U(:,idx);

        % band-basis
        JA   = U'*(JAs*U);  JB   = U'*(JBs*U);
        vA_e = U'*(vA *U);  vB_e = U'*(vB *U);

        % precompute X and denominator
        Xab(:,:,t) = 2*imag( JA .* (vB_e.') );
        Xba(:,:,t) = 2*imag( JB .* (vA_e.') );

        dE  = E(:)'-E(:);                % 4x4
        den = dE.^2 + eta^2;
        den(1:5:end) = Inf;              % diag inf
        Den(:,:,t) = den;

        Eall(:,t) = E;
    end

    % ---- pack ----
    cache.Nk    = Nk;
    cache.Norb  = Norb;
    cache.w     = w;
    cache.Xab   = Xab;
    cache.Xba   = Xba;
    cache.Den   = Den;
    cache.E     = Eall;
    cache.eta   = eta;
    cache.alpha = alpha; cache.beta = beta; cache.gamma = gamma;
    cache.hbar  = hbar;
    cache.e     = params.electronic_charge;

    % optional diagnostics at Γ
    H_G = Hk(0,0,0);
    [Vg,Dg] = eig(full(H_G));
    Eg = real(diag(Dg)); [Eg,ixg] = sort(Eg,'ascend'); Vg = Vg(:,ixg);
    Sz_exp = real(diag(Vg'*(Sg*Vg)));
    cache.sample_Gamma = [Eg, Sz_exp];
end
