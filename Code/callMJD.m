function price = callMJD(S0, K, r, q, sigma, t, gamma, delta, jumpProb)
% Call price of JD
%
% Input:
%   K: row vector
%   jumpProb: column vector
%
% Output: price(n, K)

N = numel(jumpProb) - 1;
n = (0 : N)';

muX = n * gamma;
varX = n * delta .^ 2;

eta = log(sum(jumpProb .* exp(muX + varX / 2), 1));
mu = r - q - eta / t;
muD = (mu - sigma ^ 2 / 2) * t;
varD = sigma ^ 2 * t;

muR = muD + muX;
varR = varD + varX;

sigmahat = sqrt(varR / t);
rhat = muR / t + q + sigmahat .^ 2 / 2;

price = sum(jumpProb .* exp(-(r - rhat) * t) .* callBS(S0, K, rhat, q, sigmahat, t), 1);

end

function call = callBS(S0, K, r, q, sigma, t)

d1 = (log(S0 ./ K) + (r - q + sigma .^ 2 / 2) .* t) ./ (sigma .* sqrt(t));
d2 = d1 - sigma .* sqrt(t);
nd1 = 0.5 * erfc(-d1 / sqrt(2));
nd2 = 0.5 * erfc(-d2 / sqrt(2));
call = S0 .* exp(-q .* t) .* nd1 - K .* exp(-r .* t) .* nd2;

end