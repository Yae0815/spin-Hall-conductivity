function out = spin_hall_main(params)
% σ^{s_γ}_{αβ} with Kubo (AB−BA antisymmetrized), parfor over kx
% Required params: ftn58, Ef, Nk, eta, electronic_charge, hbar, alpha, beta, gamma
% Optional: params.shift = [sx,sy,sz] in [0,1) for Monkhorst–Pack shift (grid units)

    % --- unpack ---
    %ftn   = params.ftn58;
    Ef    = params.Ef;
    Nk    = params.Nk;
    eta   = params.eta;
    e     = params.electronic_charge;
    hbar  = params.hbar;
    alpha = lower(params.alpha);
    beta  = lower(params.beta);
    gamma = lower(params.gamma);

    % --- builders / sizes ---
    %build = make_builders(ftn);
    build = make_builders_auto(params);
    Norb  = build.Norb;
    Sg    = pick_spin_op(gamma, Norb);

    % --- Γ sampling ---
    H_G = build.H(0,0,0);
    comm_SgH = Sg*H_G - H_G*Sg;
    comm_norm = norm(full(comm_SgH), 'fro');
    [Vg, Dg] = eig(full(H_G));
    Eg = real(diag(Dg));
    [Eg,ix] = sort(Eg,'ascend'); Vg = Vg(:,ix);
    Sz_exp = real(diag(Vg'*(Sg*Vg)));
    sample_Gamma = [Eg, Sz_exp];

    % --- k-mesh (support shift) ---
    base = (0:Nk-1)/Nk - 0.5;                 % periodic, no duplicate endpoints
    shift = [0 0 0];
    if isfield(params,'shift'), shift = params.shift; end   % [0,1) in grid units
    kx_list = base + shift(1)/Nk;
    ky_list = base + shift(2)/Nk;
    kz_list = base + shift(3)/Nk;
    w = 1.0 / (Nk^3);

    % --- capture handles once (broadcast into parfor) ---
    H_of_k = build.H;
    dHx    = build.dHdkx;
    dHy    = build.dHdky;
    dHz    = build.dHdkz;

    % ===== parfor over kx slices =====
    sigma_partial = zeros(Nk,1);

    parfor ixk = 1:Nk
        kx = kx_list(ixk);
        local_sum = 0.0;

        for iyk = 1:Nk
            ky = ky_list(iyk);
            for izk = 1:Nk
                kz = kz_list(izk);

                % H and velocities
                Hk  = H_of_k(kx,ky,kz);
                vA  = pick_velocity_handles(dHx,dHy,dHz, alpha, kx,ky,kz, hbar);
                vB  = pick_velocity_handles(dHx,dHy,dHz, beta,  kx,ky,kz, hbar);

                % numerical hermitization (safe)
                vA  = (vA + vA')/2;   vB  = (vB + vB')/2;
                JAs = 0.5*(Sg*vA + vA*Sg);
                JBs = 0.5*(Sg*vB + vB*Sg);
                JAs = (JAs + JAs')/2; JBs = (JBs + JBs')/2;

                % diag once
                [U,D] = eig(full(Hk));
                E = real(diag(D)); [E,idx] = sort(E,'ascend'); U = U(:,idx);

                % eigenbasis
                J_A = U'*(JAs*U);  J_B = U'*(JBs*U);
                v_A = U'*(vA *U);  v_B = U'*(vB *U);

                % occupancy T=0
                occ = (E < Ef);

                % Kubo, antisymmetrized
                Om_AB = pair(J_A, v_B, E, eta, occ);
                Om_BA = pair(J_B, v_A, E, eta, occ);
                Om = 0.5*(Om_AB - Om_BA);

                local_sum = local_sum + Om * w;
            end
        end

        sigma_partial(ixk) = local_sum;
    end

    sigma_val = (e/hbar) * sum(sigma_partial);

    % --- outputs ---
    out.sigma_sab_gamma = sigma_val;
    out.alpha = alpha; out.beta = beta; out.gamma = gamma;
    out.commutator_norm_Gamma = comm_norm;
    out.sample_bands_Gamma = sample_Gamma;
    out.note = 'lattice units with reduced k; prefactor e/ħ applied.';
end

% ===== helpers =====
function Sg = pick_spin_op(gamma, Norb)
    assert(mod(Norb,2)==0, 'Norb must be even (spin×orbital).');
    Norb_orb = Norb/2;
    switch lower(gamma)
        case 'x', s = [0 1; 1 0]/2;
        case 'y', s = [0 -1i; 1i 0]/2;
        case 'z', s = [1 0; 0 -1]/2;
        otherwise, error('gamma must be x/y/z');
    end
    Sg = kron(s, speye(Norb_orb));
end

function v = pick_velocity_handles(dHx,dHy,dHz, alpha, kx, ky, kz, hbar)
    twoPi = 2*pi;
    switch lower(alpha)
        case 'x', dH = dHx(kx,ky,kz);
        case 'y', dH = dHy(kx,ky,kz);
        case 'z', dH = dHz(kx,ky,kz);
        otherwise, error('alpha must be x/y/z');
    end
    v = (1/hbar) * (1/twoPi) * dH;
end

function Om = pair(Jmat, vmat, E, eta, occ)
    Norb = numel(E);
    Om = 0.0;
    for n = 1:Norb
        if ~occ(n), continue; end
        dE2 = (E - E(n)).^2 + eta^2;
        dE2 = max(dE2, eta^2);   % guard
        num = Jmat(n,:).*vmat(:,n).';
        num(n) = 0;
        dE2(n) = Inf;
        Om = Om + 2*imag( sum( num ./ dE2.' ) );
    end
end

function build = make_builders_auto(params)
    if isfield(params,'builder') && strcmpi(params.builder,'taguchi2020')
        assert(isfield(params,'model'), 'params.model 缺少 Taguchi 模型參數');
        build = make_builders_taguchi(params.model);   % 顯函數 DSM (Taguchi 2020)
    else
        build = make_builders(params.ftn58);           % 既有 ftn58 路線
    end
end
