function v = pick_velocity(dHx,dHy,dHz, dir, kx,ky,kz, hbar)
% shc.pick_velocity
% Unified velocity selector to avoid duplicate "speed" helpers.
% v_dir = (∂H/∂k_dir)/ħ
    switch dir
        case 'x', v = dHx(kx,ky,kz) / hbar;
        case 'y', v = dHy(kx,ky,kz) / hbar;
        case 'z', v = dHz(kx,ky,kz) / hbar;
        otherwise
            error('pick_velocity: dir must be x/y/z.');
    end
    % Hermitize to suppress tiny asymmetries
    v = (v + v')/2;
end