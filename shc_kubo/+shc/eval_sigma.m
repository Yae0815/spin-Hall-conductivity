function out = eval_sigma(cache, mu, Ef, T, method)
% shc.eval_sigma
% Evaluate σ^{s_γ}_{αβ}(μ,T) using precomputed Xab/Xba/Den/E.

    if nargin<5 || isempty(method), method = 'weighted'; end
    kB = 8.617333262e-5; % eV/K
    K   = size(cache.E, 2);
    w   = cache.w;

    sigma_sum = 0.0;

    switch lower(method)
        case 'weighted'
            for t = 1:K
                E   = cache.E(:,t);
                w_n = (T<=0) .* double(E<(mu+Ef)) + (T>0) .* (1./(1+exp((E-(mu+Ef))/(kB*T))));
                W = repmat(w_n(:), 1, numel(w_n));  % broadcast to rows
                OmAB = sum(sum( (cache.Xab(:,:,t) - cache.Xba(:,:,t)) .* (W ./ cache.Den(:,:,t)) ));
                sigma_sum = sigma_sum + real(OmAB) * w(t);
            end

        case 'bastin'
            for t = 1:K
                E   = cache.E(:,t);
                f   = (T<=0) .* double(E<(mu+Ef)) + (T>0) .* (1./(1+exp((E-(mu+Ef))/(kB*T))));
                Fnm = f(:) - f(:).';   % (f_n - f_m)
                OmAB = sum(sum( (cache.Xab(:,:,t) - cache.Xba(:,:,t)) .* (Fnm ./ cache.Den(:,:,t)) ));
                sigma_sum = sigma_sum + real(OmAB) * w(t);
            end

        otherwise
            error('Unknown method: %s', method);
    end

    out = sigma_sum;
end
