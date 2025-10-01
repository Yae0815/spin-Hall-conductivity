function build = make_builders_taguchi(params)
% builders.make_builders_taguchi
% Analytic (non-FTN58) Dirac/DSM-like 4x4 model used in Taguchi 2020 spirit.
%
% params.model requires fields: eta, txy, tz, M, beta, gamma  (all in eV)
% Example:
%   params.model = struct('eta',0.89,'txy',txy,'tz',tz,'M',M,'beta',beta,'gamma',gamma);

    mustHave = ["eta","txy","tz","M","beta","gamma"];
    assert(isfield(params,'model'), 'params.model missing');
    for f = mustHave
        assert(isfield(params.model, f), 'params.model.%s missing', f);
    end
    p = params.model;

    % Pauli
    sx = [0 1; 1 0];  sy = [0 -1i; 1i 0];  sz = [1 0; 0 -1];  s0 = eye(2);
    % Dirac-like blocks
    % Choose basis = (|↑,A>, |↓,A>, |↑,B>, |↓,B>)
    % Gamma representation (one common choice):
    %   Γ1 = sx ⊗ s0,   Γ2 = sy ⊗ s0,   Γ3 = sz ⊗ sx,   Γ4 = sz ⊗ sy,   Γ5 = sz ⊗ sz
    % You can switch to your own convention; just update H and dHdk consistently.
    Sx = kron(0.5*sx, s0);     % spin-1/2 in first factor
    Sy = kron(0.5*sy, s0);
    Sz = kron(0.5*sz, s0);

    G1 = kron(sx, s0);
    G2 = kron(sy, s0);
    G3 = kron(sz, sx);
    G4 = kron(sz, sy);
    G5 = kron(sz, sz);
    I4 = eye(4);

    % Hamiltonian components
    function Hk = H(k)
        kx = k(1); ky = k(2); kz = k(3);

        % mass-like term (txy,tz,M) （依你習慣可替換）
        m_k = p.M - p.txy*(cos(kx)+cos(ky)) - p.tz*cos(kz);

        % Dirac-like linear terms with eta; extra tilt/aniso via beta/gamma（依你習慣可替換）
        d1 = p.eta * sin(kx);         % ~ Γ1
        d2 = p.eta * sin(ky);         % ~ Γ2
        d3 = p.eta * sin(kz);         % ~ Γ3
        d4 = p.beta * sin(kz);        % ~ Γ4   (often used to split DSM to WSM/tilt)
        % particle-hole asym / additional anisotropy via gamma in-plane:
        dPHS = (p.gamma) * (2 - cos(kx) - cos(ky));  % optional PH-asym shift
        % Total
        Hk = d1*G1 + d2*G2 + d3*G3 + d4*G4 + (m_k + dPHS)*G5;
    end

    % dH/dk
    function dH = dHdk(comp, k)
        kx = k(1); ky = k(2); kz = k(3);
        switch comp
            case 'x'
                dm =  + p.txy*sin(kx);                 % d(m_k)/dkx = +txy sin kx
                dP =  + p.gamma * sin(kx);             % d(PHS)/dkx
                d1 =  + p.eta * cos(kx);               % d(d1)/dkx
                d2 =  0;
                d3 =  0;
                d4 =  0;
                dH = d1*G1 + d2*G2 + d3*G3 + d4*G4 + (dm + dP)*G5;

            case 'y'
                dm =  + p.txy*sin(ky);
                dP =  + p.gamma * sin(ky);
                d1 =  0;
                d2 =  + p.eta * cos(ky);
                d3 =  0;
                d4 =  0;
                dH = d1*G1 + d2*G2 + d3*G3 + d4*G4 + (dm + dP)*G5;

            case 'z'
                dm =  + p.tz*sin(kz);
                dP =  0;
                d1 =  0;
                d2 =  0;
                d3 =  + p.eta * cos(kz);
                d4 =  + p.beta * cos(kz);
                dH = d1*G1 + d2*G2 + d3*G3 + d4*G4 + (dm + dP)*G5;

            otherwise
                error('dHdk: unknown component %s', comp);
        end
    end

    function V = vel(comp, k, hbar)
        if nargin<3 || isempty(hbar), hbar = 1; end
        V = (1/hbar) * dHdk(comp, k);
    end

    function S = spin_op(which)
        switch lower(which)
            case {'x','sx'}, S = Sx;
            case {'y','sy'}, S = Sy;
            case {'z','sz'}, S = Sz;
            otherwise, error('spin_op: unknown %s', which);
        end
    end

    build = struct();
    build.Norb = 4;
    build.H    = @H;
    build.dHdk = @dHdk;
    build.vel  = @vel;
    build.spin = @spin_op;
end
