function out = eval_sigma(cache, mu, Ef, T, method)
% shc.eval_sigma
% Evaluate σ^{s_γ}_{αβ}(μ,T) using precomputed Xab/Xba/Den/E.
% method: 'weighted' (Fermi-sea w_n) or 'bastin' (f_n - f_m)

    % ---- defaults (in order) ----
    if nargin < 3 || isempty(Ef),     Ef = 0;          end
    if nargin < 4 || isempty(T),      T  = 0;          end
    if nargin < 5 || isempty(method), method = 'weighted'; end

    kB = 8.617333262e-5; % eV/K
    [Norb, K] = size(cache.E);

    w   = cache.w;
    sigma_sum = 0.0;

    switch lower(method)
        case 'weighted'
            if T <= 0
                % T=0 → Fermi-sea projector
                for t = 1:K
                    E   = cache.E(:,t);
                    w_n = double(E < (mu+Ef));            % Norb×1
                    W   = w_n(:) * ones(1, Norb);         % Norb×Norb
                    OmAB = sum(sum( (cache.Xab(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                    OmBA = sum(sum( (cache.Xba(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                    sigma_sum = sigma_sum + 0.5*(OmAB - OmBA) * w;
                end
            else
                % finite T → Fermi-Dirac weights
                invkBT = 1/(kB*T);
                for t = 1:K
                    E   = cache.E(:,t);
                    fFD = 1 ./ (1 + exp((E - (mu+Ef)) * invkBT));  % Norb×1
                    W   = fFD(:) * ones(1, Norb);                  % Norb×Norb
                    OmAB = sum(sum( (cache.Xab(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                    OmBA = sum(sum( (cache.Xba(:,:,t) ./ cache.Den(:,:,t)) .* W ));
                    sigma_sum = sigma_sum + 0.5*(OmAB - OmBA) * w;
                end
            end

        case 'bastin'
            for t = 1:K
                E = cache.E(:,t);
                if T <= 0
                    f = double(E < (mu+Ef));
                else
                    f = 1 ./ (1 + exp((E - (mu+Ef)) / (kB*T)));
                end
                F = f(:) - f(:)';          % Norb×Norb
                F(1:Norb+1:end) = 0;       % 清掉對角 (避免 0/Inf 的不必要操作)
                OmAB = sum(sum( (cache.Xab(:,:,t) ./ cache.Den(:,:,t)) .* F ));
                OmBA = sum(sum( (cache.Xba(:,:,t) ./ cache.Den(:,:,t)) .* F ));
                sigma_sum = sigma_sum + 0.5*(OmAB - OmBA) * w;
            end

        otherwise
            error('method must be ''weighted'' or ''bastin''.');
    end

    out.sigma = (cache.e/cache.hbar) * sigma_sum;
end
