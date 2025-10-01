function build = make_builders_taguchi(model)
% make_builders_taguchi
%   Taguchi et al. (2020) Topological Dirac semimetal lattice model:
%   H(k) = a1 σx sz + a2 σy s0 + a3 σz s0 + a4 σx sx + a5 σx sy
%   with
%     a1 =  η sin kx
%     a2 = -η sin ky
%     a3 =  M - txy (cos kx + cos ky) - tz cos kz
%     a4 = (β+γ) sin kz (cos ky - cos kx)
%     a5 = -(β-γ) sin kz sin ky sin kx
%   k 以**物理座標**（-π..π）定義；外部呼叫仍使用 reduced k，所以此函式內部會做 k = 2π*k_red。
%
%   Basis: (s↑, p↑, s↓, p↓)  ⇒  ordering = spin ⊗ orbital
%   σi: orbital Pauli; si: spin Pauli
%
%   Inputs (model): struct with fields
%       eta, txy, tz, beta, gamma, M
%
%   Output (build): struct with handles
%       H, dHdkx, dHdky, dHdkz, Norb

    % ---- unpack parameters ----
    eta   = model.eta;
    txy   = model.txy;
    tz    = model.tz;
    beta  = model.beta;
    gamma = model.gamma;
    M0    = model.M;

    Norb = 4;

    % ---- Pauli (spin ⊗ orbital) ----
    s0 = eye(2); sx = [0 1;1 0]; sy = [0 -1i; 1i 0]; sz = [1 0;0 -1];
    soI = @(A) kron(A, eye(2));   % spin ⊗ I_orb
    Is  = @(A) kron(eye(2), A);   % I_spin ⊗ orbital
    % Building blocks appearing in H:
    X_sy0 = Is([0 -1i; 1i 0]);          % σy s0  = I_s ⊗ σy
    X_sz0 = Is([1 0; 0 -1]);            % σz s0  = I_s ⊗ σz
    X_sxsz= kron(sz, [0 1; 1 0]);       % σx sz  = sz ⊗ σx
    X_sxsy= kron(sy, [0 1; 1 0]);       % σx sy  = sy ⊗ σx
    X_sxsx= kron(sx, [0 1; 1 0]);       % σx sx  = sx ⊗ σx

    twoPi = 2*pi;

    % ---- scalar ai(k) + their partials (k here is physical) ----
    function [a1,a2,a3,a4,a5] = a_all(kx,ky,kz)
        a1 =  eta * sin(kx);
        a2 = -eta * sin(ky);
        a3 =  M0 - txy*(cos(kx)+cos(ky)) - tz*cos(kz);
        a4 =  (beta+gamma) * sin(kz) * (cos(ky) - cos(kx));
        a5 = -(beta-gamma) * sin(kz) * sin(ky) * sin(kx);
    end

    function [d1x,d2x,d3x,d4x,d5x] = dax(kx,ky,kz)  % ∂/∂kx
        d1x =  eta * cos(kx);
        d2x =  0;
        d3x =  txy * sin(kx);
        d4x =  (beta+gamma) * sin(kz) * (+sin(kx));           % d/dkx[-cos(kx)] = +sin(kx)
        d5x = -(beta-gamma) * sin(kz) * sin(ky) * cos(kx);
    end

    function [d1y,d2y,d3y,d4y,d5y] = day(kx,ky,kz)  % ∂/∂ky
        d1y =  0;
        d2y = -eta * cos(ky);
        d3y =  txy * sin(ky);
        d4y =  (beta+gamma) * sin(kz) * (-sin(ky));           % d/dky[+cos(ky)] = -sin(ky)
        d5y = -(beta-gamma) * sin(kz) * cos(ky) * sin(kx);
    end

    function [d1z,d2z,d3z,d4z,d5z] = daz(kx,ky,kz)  % ∂/∂kz
        d1z = 0;
        d2z = 0;
        d3z = tz * sin(kz);
        d4z = (beta+gamma) * cos(kz) * (cos(ky) - cos(kx));
        d5z = -(beta-gamma) * cos(kz) * sin(ky) * sin(kx);
    end

    % ---- H(k_red) & dH/dk_red (include 2π chain-rule) ----
    function Hk = H_of_k(kx_red,ky_red,kz_red)
        kx = twoPi*kx_red; ky = twoPi*ky_red; kz = twoPi*kz_red;
        [a1,a2,a3,a4,a5] = a_all(kx,ky,kz);
        Hk = a1*X_sxsz + a2*X_sy0 + a3*X_sz0 + a4*X_sxsx + a5*X_sxsy;
        Hk = (Hk + Hk')/2;             % 保險性 Hermitize
        Hk = sparse(Hk);
    end

    function dHx = dHdkx(kx_red,ky_red,kz_red)
        kx = twoPi*kx_red; ky = twoPi*ky_red; kz = twoPi*kz_red;
        [d1,d2,d3,d4,d5] = dax(kx,ky,kz);    % ∂/∂kx (physical)
        dH = d1*X_sxsz + d2*X_sy0 + d3*X_sz0 + d4*X_sxsx + d5*X_sxsy;
        dHx = twoPi * sparse( (dH + dH')/2 ); % ∂/∂k_red = 2π * ∂/∂kx
    end

    function dHy = dHdky(kx_red,ky_red,kz_red)
        kx = twoPi*kx_red; ky = twoPi*ky_red; kz = twoPi*kz_red;
        [d1,d2,d3,d4,d5] = day(kx,ky,kz);
        dH = d1*X_sxsz + d2*X_sy0 + d3*X_sz0 + d4*X_sxsx + d5*X_sxsy;
        dHy = twoPi * sparse( (dH + dH')/2 );
    end

    function dHz = dHdkz(kx_red,ky_red,kz_red)
        kx = twoPi*kx_red; ky = twoPi*ky_red; kz = twoPi*kz_red;
        [d1,d2,d3,d4,d5] = daz(kx,ky,kz);
        dH = d1*X_sxsz + d2*X_sy0 + d3*X_sz0 + d4*X_sxsx + d5*X_sxsy;
        dHz = twoPi * sparse( (dH + dH')/2 );
    end

    build.H     = @H_of_k;
    build.dHdkx = @dHdkx;
    build.dHdky = @dHdky;
    build.dHdkz = @dHdkz;
    build.Norb  = Norb;
end
