function out = eval_sigma(cache, mu, Ef, T, method)
% shc.eval_sigma
% Evaluate σ^{s_γ}_{αβ}(μ,T) using precomputed Xab/Xba/Den/E.
% method: 'weighted' (Fermi-sea w_n) or 'bastin' (f_n - f_m)

    if nargin<4 || isempty(method), method = 'weighted'; end
    kB = 8.617333262e-5; % eV/K
    K   = size(cache.E, 2);
    w   = cache.w;

    sigma_sum = 0.0;

    switch lower(method)
        case 'weighted'
            for t = 1:K
                E   = cache.E(:,t);
                w_n = (T<=0) .* double(E<(mu+Ef)) + (T>0) .* (1./(1+exp((E-(mu+Ef))/(kB*T))));
                % broadcast w_n to rows
                W = repmat(w_n(:), 1, 4);  % 4x4
                OmAB = sum(sum( (cache.Xab(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                OmBA = sum(sum( (cache.Xba(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                sigma_sum = sigma_sum + 0.5*(OmAB - OmBA) * w;
            end

        case 'bastin'
            for t = 1:K
                E = cache.E(:,t);
                if T<=0
                    f = double(E<(mu+Ef));
                else
                    f = 1./(1+exp((E-(mu+Ef))/(kB*T)));
                end
                F = f(:) - f(:)';          % 4x4
                F(1:5:end) = 0;            % remove diag
                OmAB = sum(sum( (cache.Xab(:,:,t) ./ cache.Den(:,:,t)) .* F ));
                OmBA = sum(sum( (cache.Xba(:,:,t) ./ cache.Den(:,:,t)) .* F ));
                sigma_sum = sigma_sum + 0.5*(OmAB - OmBA) * w;
            end

        otherwise
            error('method must be ''weighted'' or ''bastin''.');
    end

    out.sigma = (cache.e/cache.hbar) * sigma_sum;
end
