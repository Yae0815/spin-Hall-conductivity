function build = shc_builders_dsm_taguchi(params)
% returns function handles for H(k), dHdk, Sg
    eta  = params.eta_band;   % <-- 用於 a1,a2 的 η（與 Kubo 展寬不同名）
    txy  = params.txy;
    tz   = params.tz;
    beta = params.beta;
    gamma= params.gamma;
    M    = params.M;

    % Pauli
    sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1]; s0 = eye(2);
    % tensor helpers
    kron2 = @(A,B) kron(A,B);

    function H = Hk(kx,ky,kz)
        a1 = eta*sin(kx);
        a2 = -eta*sin(ky);
        a3 = M - txy*(cos(kx)+cos(ky)) - tz*cos(kz);
        a4 = (beta+gamma)*sin(kz)*(cos(ky)-cos(kx));
        a5 = -(beta-gamma)*sin(kz)*sin(ky)*sin(kx);
        H  = a1*kron2(sx,sz) + a2*kron2(sy,s0) + a3*kron2(sz,s0) ...
           + a4*kron2(sx,sx) + a5*kron2(sx,sy);
    end

    function d = dHdkx(kx,ky,kz)
        d =  eta*cos(kx)*kron2(sx,sz) ...
           + txy*sin(kx)*kron2(sz,s0) ...
           + (beta+gamma)*sin(kz)*sin(kx)*kron2(sx,-sx) ...
           - (beta-gamma)*sin(kz)*sin(ky)*cos(kx)*kron2(sx,sy);
    end
    function d = dHdky(kx,ky,kz)
        d = -eta*cos(ky)*kron2(sy,s0) ...
           + txy*sin(ky)*kron2(sz,s0) ...
           + (beta+gamma)*sin(kz)*(-sin(ky))*kron2(sx,sx) ...
           - (beta-gamma)*sin(kz)*cos(ky)*sin(kx)*kron2(sx,sy);
    end
    function d = dHdkz(kx,ky,kz)
        d =  tz*sin(kz)*kron2(sz,s0) ...
           + (beta+gamma)*cos(kz)*(cos(ky)-cos(kx))*kron2(sx,sx) ...
           - (beta-gamma)*cos(kz)*sin(ky)*sin(kx)*kron2(sx,sy);
    end

    function S = Sg_sel(gamma_char)
        switch lower(gamma_char)
          case 'z', S = 0.5*[1 0;0 -1]; % ℏ/2 absorbed later if you keep ħ=1
                 S = kron2(s0, [1 0;0 -1]); % spin on 2nd space
          case 'x', S = kron2(s0, sx);
          case 'y', S = kron2(s0, sy);
        end
    end

    build.H      = @Hk;
    build.dHdkx  = @dHdkx;
    build.dHdky  = @dHdky;
    build.dHdkz  = @dHdkz;
    build.Sx     = Sg_sel('x'); build.Sy = Sg_sel('y'); build.Sz = Sg_sel('z');
    build.Norb   = 4;
end
