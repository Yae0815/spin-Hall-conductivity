function S = util_pick_spin_op(build, gamma)
% shc.util_pick_spin_op
    if isa(build.spin, 'function_handle')
        S = build.spin(gamma);
    else
        error('build.spin must be a function handle returning Sx/Sy/Sz');
    end
end

