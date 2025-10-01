function build = make_builders(params)
% shc.make_builders
% Return a struct "build" providing H(k), dHdk, velocity, spin op, etc.
% Compatible with both analytic (non-FTN58) and future FTN58 backends.

    assert(isfield(params,'builder'), 'params.builder is required');

    if isa(params.builder, 'function_handle')
        build = params.builder(params);
        return;
    end

    switch lower(string(params.builder))
        case "taguchi2020"
            build = builders.make_builders_taguchi(params);

        case "ftn58"
            % Optional: only if/when you want to hook FTN58 back
            build = builders.make_builders_ftn58(params);

        otherwise
            error('Unknown builder: %s', string(params.builder));
    end
end
