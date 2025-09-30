function build = shc_builders_dsm_taguchi(params)
% shc.shc_builders_dsm_taguchi
% Returns function handles for H(k), its derivatives, and spin operators.
% Base conventions:
%   - Norb = 4
%   - Basis = (spin ⊗ orbital) 2x2
%   - Spin operators act on the FIRST space (kron(σ, I))
%
% Model parameters (use *_cpl to avoid name collision with tensor indices):
%   params.eta_band, params.txy, params.tz, params.beta_cpl, params.gamma_cpl, params.M

    eta   = params.eta_band;
    txy   = params.txy;
    tz    = params.tz;
    betaC = params.beta_cpl;     % coupling "beta"
    gammaC= params.gamma_cpl;    % coupling "gamma"
    M     = params.M;

    % Pauli
    sx = [0 1; 1 0]; 
    sy = [0 -1i; 1i 0]; 
    sz = [1 0; 0 -1]; 
    s0 = eye(2);

    kron2 = @(A,B) kron(A,B);

    % ----- H(k) -----
    function H = Hk(kx,ky,kz)
        a1 =  eta * sin(kx);
        a2 = -eta * sin(ky);
        a3 =  M - txy*(cos(kx)+cos(ky)) - tz*cos(kz);
        a4 = (betaC+gammaC)*sin(kz)*(cos(ky)-cos(kx));
        a5 = -(betaC-gammaC)*sin(kz)*sin(ky)*sin(kx);

        % spin in first space, orbital in second
        H  = a1*kron2(sx,sz) + a2*kron2(sy,s0) + a3*kron2(sz,s0) ...
           + a4*kron2(sx,sx) + a5*kron2(sx,sy);
    end

    % ----- dH/dk -----
    function d = dHdkx(kx,ky,kz)
        d =  eta*cos(kx)*kron2(sx,sz) ...
           + txy*sin(kx)*kron2(sz,s0) ...
           + (betaC+gammaC)*sin(kz)*sin(kx)*kron2(sx,-sx) ...
           - (betaC-gammaC)*sin(kz)*sin(ky)*cos(kx)*kron2(sx,sy);
    end

    function d = dHdky(kx,ky,kz)
        d = -eta*cos(ky)*kron2(sy,s0) ...
           + txy*sin(ky)*kron2(sz,s0) ...
           + (betaC+gammaC)*sin(kz)*(-sin(ky))*kron2(sx,sx) ...
           - (betaC-gammaC)*sin(kz)*cos(ky)*sin(kx)*kron2(sx,sy);
    end

    function d = dHdkz(kx,ky,kz)
        d =  tz*sin(kz)*kron2(sz,s0) ...
           + (betaC+gammaC)*cos(kz)*(cos(ky)-cos(kx))*kron2(sx,sx) ...
           - (betaC-gammaC)*cos(kz)*sin(ky)*sin(kx)*kron2(sx,sy);
    end

    % ----- Spin operators: spin ⊗ I (include 1/2) -----
    Sx = 0.5*kron2(sx, s0);
    Sy = 0.5*kron2(sy, s0);
    Sz = 0.5*kron2(sz, s0);

    build.H      = @Hk;
    build.dHdkx  = @dHdkx;
    build.dHdky  = @dHdky;
    build.dHdkz  = @dHdkz;

    build.Sx     = Sx; 
    build.Sy     = Sy; 
    build.Sz     = Sz;
    build.S.x    = Sx; 
    build.S.y    = Sy; 
    build.S.z    = Sz;

    build.Norb   = 4;
end
