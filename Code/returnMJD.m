function [f, F] = returnMJD(y, r, q, sigma, t, gamma, delta, jumpProb)
% Return PDF and CDF of JD
%
% Input:
%   y: row vector
%   jumpProb: column vector
%
% Output:
%   fR: column vector

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

f = sum(jumpProb .* exp(-(y - muR) .^2 ./ (2 * varR)) ./ sqrt(2 * pi .* varR), 1)';
F = sum(jumpProb * 0.5 .* erfc(-(y - muR) ./ sqrt(2 * varR)), 1)';

end
