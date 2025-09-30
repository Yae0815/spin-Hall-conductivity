function cache = precompute_kgrid(params)
% shc.precompute_kgrid
% One-time preprocessing for Spin Hall (Kubo AB-BA) with N=4.
% Stores, for each k, E, Xab, Xba, and denom matrix.
% Changing eta/alpha/beta/gamma/Nk/shift/ftn requires re-precompute.

    % ---- unpack ----
    %Sanity Check
    %ftn   = params.ftn58;
    Nk    = params.Nk;
    eta   = params.eta;
    hbar  = params.hbar;
    alpha = lower(params.alpha);
    beta  = lower(params.beta);
    gamma = lower(params.gamma);

    % ---- builders ----
    %Sanity Check
    %build = shc.make_builders(ftn);
    build = shc.shc_builders_dsm_taguchi(params);
    Norb  = build.Norb;
    assert(Norb==4, 'This precompute is optimized for Norb=4.');
    Sg    = loc_pick_spin_op(gamma, Norb);

    % ---- k grid ----
    base = (0:Nk-1)/Nk - 0.5;
    shift = [0 0 0];
    if isfield(params,'shift'), shift = params.shift; end
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

    % ---- flatten k loop for better parfor balance ----
    % map linear idx -> (ix,iy,iz)
    [IX,IY,IZ] = ndgrid(1:Nk, 1:Nk, 1:Nk);
    ix = IX(:); iy = IY(:); iz = IZ(:);

    parfor t = 1:K
        kx = kx_list(ix(t)); ky = ky_list(iy(t)); kz = kz_list(iz(t));

        % H and velocities
        H   = Hk(kx,ky,kz);
        vA  = loc_pick_velocity(dHx,dHy,dHz, alpha, kx,ky,kz, hbar);
        vB  = loc_pick_velocity(dHx,dHy,dHz, beta,  kx,ky,kz, hbar);

        % hermitize
        vA  = (vA + vA')/2; vB = (vB + vB')/2;
        JAs = 0.5*(Sg*vA + vA*Sg);
        JBs = 0.5*(Sg*vB + vB*Sg);
        JAs = (JAs + JAs')/2; JBs = (JBs + JBs')/2;

        % diag & band-basis
        [U,D] = eig(full(H));
        E = real(diag(D)); [E,idx] = sort(E,'ascend'); U = U(:,idx);

        JA = U'*(JAs*U);  JB = U'*(JBs*U);
        vA_e = U'*(vA *U); vB_e = U'*(vB *U);

        % precompute X and denominator
        Xab(:,:,t) = 2*imag( JA .* (vB_e.') );
        Xba(:,:,t) = 2*imag( JB .* (vA_e.') );

        dE = E(:)'-E(:);                       % 4x4
        den = dE.^2 + eta^2;
        den(1:5:end) = Inf;                    % diag inf
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
    H_G = build.H(0,0,0);
    [Vg,Dg] = eig(full(H_G));
    Eg = real(diag(Dg)); [Eg,ixg] = sort(Eg,'ascend'); Vg = Vg(:,ixg);
    Sz_exp = real(diag(Vg'*(Sg*Vg)));
    cache.sample_Gamma = [Eg, Sz_exp];
end

% ==== local helpers =======================================================
function Sg = loc_pick_spin_op(gamma, Norb)
    assert(mod(Norb,2)==0);
    Norb_orb = Norb/2;
    switch lower(gamma)
        case 'x', s = [0 1; 1 0]/2;
        case 'y', s = [0 -1i; 1i 0]/2;
        case 'z', s = [1 0; 0 -1]/2;
        otherwise, error('gamma must be x/y/z');
    end
    Sg = kron(s, speye(Norb_orb));
end

function v = loc_pick_velocity(dHx,dHy,dHz, alpha, kx, ky, kz, hbar)
    twoPi = 2*pi;
    switch lower(alpha)
        case 'x', dH = dHx(kx,ky,kz);
        case 'y', dH = dHy(kx,ky,kz);
        case 'z', dH = dHz(kx,ky,kz);
        otherwise, error('alpha must be x/y/z');
    end
    v = (1/hbar) * (1/twoPi) * dH;
end
